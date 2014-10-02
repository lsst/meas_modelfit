// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#include <cstdlib>

#include "boost/filesystem/path.hpp"
#include "boost/make_shared.hpp"

#include "ndarray/eigen.h"

#include "lsst/afw/detection/FootprintSet.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/shapelet/FunctorKeys.h"
#include "lsst/meas/multifit/TruncatedGaussian.h"
#include "lsst/meas/multifit/MultiModel.h"
#include "lsst/meas/multifit/CModel.h"
#include "lsst/meas/base/Results.h"

namespace lsst { namespace meas { namespace multifit {


//-------------------- Utility code -------------------------------------------------------------------------

namespace {

PTR(afw::detection::Footprint) _mergeFootprints(
    afw::detection::Footprint const & a,
    afw::detection::Footprint const & b
) {
    // Yuck: we have no routine that merges Footprints, so as a workaround we make a Mask from
    // the footprints, then run detection on the Mask to make a FootprintSet, then extract the
    // first Footprint from it, after checking that it's the only one.
    afw::geom::Box2I bbox(a.getBBox());
    bbox.include(b.getBBox());
    afw::image::Mask<> mask(bbox);
    afw::detection::setMaskFromFootprint(&mask, a, afw::image::MaskPixel(0x1));
    afw::detection::setMaskFromFootprint(&mask, b, afw::image::MaskPixel(0x1));
    afw::detection::FootprintSet fpSet(
        mask,
        afw::detection::Threshold(0x1, afw::detection::Threshold::BITMASK),
        1 // npixMin
    );
    if (fpSet.getFootprints()->size() > 1u) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeError,
            (boost::format("Footprints to be merged do not overlap: %s vs. %s")
             % a.getBBox() % b.getBBox()).str()
        );
    }
    return fpSet.getFootprints()->front();
}

Pixel computeFluxInFootprint(
    afw::image::Image<Pixel> const & image,
    afw::detection::Footprint const & footprint
) {
    return flattenArray(footprint, image.getArray(), image.getXY0()).asEigen().sum();
}

} // anonymous

//-------------------- Control Objects ----------------------------------------------------------------------

PTR(Model) CModelStageControl::getModel() const {
    return Model::make(getProfile().getBasis(nComponents, maxRadius), Model::FIXED_CENTER);
}

PTR(Prior) CModelStageControl::getPrior() const {
    if (priorSource == "NONE") {
        return PTR(Prior)();
    } else if (priorSource == "FILE") {
        char const * pkgDir = std::getenv("MEAS_MULTIFIT_DIR");
        if (!pkgDir) {
            throw LSST_EXCEPT(
                meas::base::FatalAlgorithmError,
                "MEAS_MULTIFIT_DIR environment variable not defined; cannot find persisted Priors"
            );
        }
        boost::filesystem::path priorPath
            = boost::filesystem::path(pkgDir)
            / boost::filesystem::path("data")
            / boost::filesystem::path(priorName + ".fits");
        PTR(Mixture) mixture = Mixture::readFits(priorPath.string());
        return boost::make_shared<MixturePrior>(mixture, "single-ellipse");
    } else if (priorSource == "CONFIG") {
        return boost::make_shared<SoftenedLinearPrior>(priorConfig);
    } else {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "priorSource must be one of 'NONE', 'FILE', or 'CONFIG'"
        );
    }
}

// ------------------- Result Objects -----------------------------------------------------------------------

CModelStageResult::CModelStageResult() :
    flux(std::numeric_limits<Scalar>::quiet_NaN()),
    fluxSigma(std::numeric_limits<Scalar>::quiet_NaN()),
    objective(std::numeric_limits<Scalar>::quiet_NaN()),
    ellipse(std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::quiet_NaN(),
            std::numeric_limits<Scalar>::quiet_NaN(), false)
{
    flags[FAILED] = true;
}

CModelResult::CModelResult() :
    flux(std::numeric_limits<Scalar>::quiet_NaN()),
    fluxSigma(std::numeric_limits<Scalar>::quiet_NaN()),
    fracDev(std::numeric_limits<Scalar>::quiet_NaN()),
    objective(std::numeric_limits<Scalar>::quiet_NaN())
{
    flags[FAILED] = true;
}


// ------------------- Key Objects for transferring to/from afw::table Records ------------------------------

namespace {

// Keys for a single stage (one of 'initial', 'exp', 'dev'
struct CModelStageKeys {

    // this constructor is used to allocate output fields in both forced and non-forced mode
    CModelStageKeys(
        Model const & model,
        afw::table::Schema & schema,
        std::string const & prefix,
        std::string const & stage,
        bool isForced,
        CModelStageControl const & ctrl
    ) :
        flux(
            schema.addField<meas::base::Flux>(
                schema.join(prefix, "flux"),
                "flux from the " + stage + " fit",
                "dn"
            )
        ),
        fluxSigma(
            schema.addField<meas::base::FluxErrElement>(
                schema.join(prefix, "fluxSigma"),
                "flux uncertainty from the " + stage + " fit",
                "dn"
            )
        ),
        fluxFlag(
            schema.addField<afw::table::Flag>(
                schema.join(prefix, "flag"),
                "flag set when the flux for the " + stage + " flux failed"
            )
        )
    {
        if (!isForced) {
            ellipse = afw::table::QuadrupoleKey::addFields(
                schema,
                schema.join(prefix, "ellipse"),
                "half-light ellipse of the " + stage + " fit",
                "pixels^2"
            );
            objective = schema.addField<Scalar>(
                schema.join(prefix, "objective"),
                "-ln(likelihood*prior) at best-fit point for the " + stage + " fit"
            );
            nonlinear = afw::table::ArrayKey<double>::addFields(
                schema,
                schema.join(prefix, "nonlinear"),
                "nonlinear parameters for the " + stage + " fit", "",
                model.getNonlinearDim()
            );
            fixed = afw::table::ArrayKey<double>::addFields(
                schema,
                schema.join(prefix, "fixed"),
                "fixed parameters for the " + stage + " fit", "",
                model.getFixedDim()
            );
            flags[CModelStageResult::TR_SMALL] = schema.addField<afw::table::Flag>(
                schema.join(prefix, "flag", "trSmall"),
                "the optimizer converged because the trust radius became too small; this is a less-secure "
                "result than when the gradient is below the threshold, but usually not a problem"
            );
            flags[CModelStageResult::MAX_ITERATIONS] = schema.addField<afw::table::Flag>(
                schema.join(prefix, "flag", "maxIter"),
                "the optimizer hit the maximum number of iterations and did not converge"
            );
            flags[CModelStageResult::NUMERIC_ERROR] = schema.addField<afw::table::Flag>(
                schema.join(prefix, "flag", "numericError"),
                "numerical underflow or overflow in model evaluation; usually this means the prior was "
                "insufficient to regularize the fit"
            );
            if (ctrl.doRecordHistory) {
                nIter = schema.addField<int>(
                    schema.join(prefix, "nIter"),
                    "Number of total iterations in stage"
                );
            }
            if (ctrl.doRecordTime) {
                time = schema.addField<Scalar>(
                    schema.join(prefix, "time"),
                    "Time spent in stage", "seconds"
                );
            }
        }
        flags[CModelStageResult::FAILED] = fluxFlag; // these flags refer to the same underlying field
    }

    // this constructor is used to get needed keys from the reference schema in forced mode
    CModelStageKeys(
        Model const & model,
        afw::table::Schema const & schema,
        std::string const & prefix
    ) :
        flux(schema[prefix]["flux"]),
        fluxSigma(schema[prefix]["fluxSigma"]),
        fluxFlag(schema[prefix]["flag"]),
        nonlinear(schema[prefix]["nonlinear"]),
        fixed(schema[prefix]["fixed"])
    {
        flags[CModelStageResult::FAILED] = fluxFlag; // these flags refer to the same underlying field
        LSST_THROW_IF_NE(
            model.getNonlinearDim(), nonlinear.getSize(),
            meas::base::FatalAlgorithmError,
            "Configured model nonlinear dimension (%d) does not match reference schema (%d)"
        );
        LSST_THROW_IF_NE(
            model.getFixedDim(), fixed.getSize(),
            meas::base::FatalAlgorithmError,
            "Configured model fixed dimension (%d) does not match reference schema (%d)"
        );
    }

    void copyResultToRecord(CModelStageResult const & result, afw::table::BaseRecord & record) {
        record.set(flux, result.flux);
        record.set(fluxSigma, result.fluxSigma);
        record.set(fluxFlag, result.getFlag(CModelStageResult::FAILED));
        if (objective.isValid()) {
            record.set(objective, result.objective);
        }
        if (ellipse.isValid()) {
            record.set(ellipse, result.ellipse);
        }
        if (nonlinear.isValid() && !result.nonlinear.isEmpty()) {
            record.set(nonlinear, result.nonlinear);
        }
        if (fixed.isValid() && !result.fixed.isEmpty()) {
            record.set(fixed, result.fixed);
        }
        if (nIter.isValid()) {
            record.set(nIter, result.history.size());
        }
        if (time.isValid()) {
            record.set(time, result.history.size());
        }
        for (int b = 0; b < CModelStageResult::N_FLAGS; ++b) {
            if (flags[b].isValid()) {
                record.set(flags[b], result.flags[b]);
            }
        }
    }

    CModelStageResult copyRecordToResult(afw::table::BaseRecord const & record) const {
        // this is only used when reading reference records, so we only transfer the fields we need for that
        CModelStageResult result;
        result.setFlag(CModelStageResult::FAILED, record.get(flags[CModelStageResult::FAILED]));
        result.nonlinear = record.get(nonlinear);
        result.fixed = record.get(fixed);
        return result;
    }

    afw::table::Key<meas::base::Flux> flux;
    afw::table::Key<meas::base::FluxErrElement> fluxSigma;
    afw::table::Key<afw::table::Flag> fluxFlag;
    afw::table::QuadrupoleKey ellipse;
    afw::table::Key<Scalar> objective;
    afw::table::Key<afw::table::Flag> flags[CModelStageResult::N_FLAGS];
    afw::table::ArrayKey<Scalar> nonlinear;
    afw::table::ArrayKey<Scalar> fixed;
    afw::table::Key<Scalar> time;
    afw::table::Key<int> nIter;
};

// Master Keys object for CModel; holds keys that aren't specific to one nonlinear stage
struct CModelKeys {

    // this constructor is used to allocate output fields in both forced and non-forced mode
    CModelKeys(
        Model const & initialModel, Model const & expModel, Model const & devModel,
        afw::table::Schema & schema,
        std::string const & prefix,
        bool isForced,
        CModelControl const & ctrl
    ) :
        initial(initialModel, schema, schema.join(prefix, "initial"), "initial",
                isForced, ctrl.initial),
        exp(expModel, schema, schema.join(prefix, "exp"), "exponential", isForced, ctrl.exp),
        dev(devModel, schema, schema.join(prefix, "dev"), "de Vaucouleur", isForced, ctrl.dev),
        // Unlike all the other keys, we expect the psf keys to already be present in the schema,
        // and we just retrieve them, because they're created and filled by the ShapeletPsfApprox plugin.
        psf(schema["multifit"]["ShapeletPsfApprox"][ctrl.psfName]),
        flux(
            schema.addField<meas::base::Flux>(
                schema.join(prefix, "flux"),
                "flux from the final cmodel fit"
            )
        ),
        fluxSigma(
            schema.addField<meas::base::FluxErrElement>(
                schema.join(prefix, "fluxSigma"),
                "flux uncertainty from the final cmodel fit"
            )
        ),
        fluxFlag(
            schema.addField<afw::table::Flag>(
                schema.join(prefix, "flag"),
                "flag set if the final cmodel fit (or any previous fit) failed"
            )
        ),
        fracDev(
            schema.addField<Scalar>(
                schema.join(prefix, "fracDev"),
                "fraction of flux in de Vaucouleur component"
            )
        ),
        objective(
            schema.addField<Scalar>(
                schema.join(prefix, "objective"),
                "-ln(likelihood) (chi^2) in cmodel fit"
            )
        )
    {
        flags[CModelResult::FAILED] = fluxFlag; // these keys refer to the same underlying field
        flags[CModelResult::MAX_AREA] = schema.addField<afw::table::Flag>(
            schema.join(prefix, "flag", "maxArea"),
            "number of pixels in fit region exceeded the region.maxArea value (usually due to bad moments)"
        );
        flags[CModelResult::MAX_BAD_PIXEL_FRACTION] = schema.addField<afw::table::Flag>(
            schema.join(prefix, "flag", "maxBadPixelFraction"),
            "the fraction of bad/clipped pixels in the fit region exceeded region.maxBadPixelFraction"
        );
        flags[CModelResult::NO_SHAPE] = schema.addField<afw::table::Flag>(
            schema.join(prefix, "flag", "noShape"),
            "the shape slot needed to initialize the parameters failed or was not defined"
        );
        flags[CModelResult::NO_SHAPELET_PSF] = schema.addField<afw::table::Flag>(
            schema.join(prefix, "flag", "noShapeletPsf"),
            "the multishapelet fit to the PSF model did not succeed"
        );
    }

    // this constructor is used to get needed keys from the reference schema in forced mode
    CModelKeys(
        Model const & initialModel, Model const & expModel, Model const & devModel,
        afw::table::Schema const & schema,
        std::string const & prefix
    ) :
        initial(initialModel, schema, schema.join(prefix, "initial")),
        exp(expModel, schema, schema.join(prefix, "exp")),
        dev(devModel, schema, schema.join(prefix, "dev"))
    {}

    void copyResultToRecord(CModelResult const & result, afw::table::BaseRecord & record) {
        initial.copyResultToRecord(result.initial, record);
        exp.copyResultToRecord(result.exp, record);
        dev.copyResultToRecord(result.dev, record);
        record.set(flux, result.flux);
        record.set(fluxSigma, result.fluxSigma);
        record.set(fracDev, result.fracDev);
        record.set(objective, result.objective);
        for (int b = 0; b < CModelResult::N_FLAGS; ++b) {
            if (flags[b].isValid()) {
                record.set(flags[b], result.flags[b]);
            }
        }
    }

    CModelResult copyRecordToResult(afw::table::BaseRecord const & record) const {
        // this is only used when reading reference records, so we only transfer the fields we need for that
        CModelResult result;
        result.initial = initial.copyRecordToResult(record);
        result.exp = exp.copyRecordToResult(record);
        result.dev = dev.copyRecordToResult(record);
        return result;
    }

    CModelStageKeys initial;
    CModelStageKeys exp;
    CModelStageKeys dev;
    shapelet::MultiShapeletFunctionKey psf;
    afw::table::Key<meas::base::Flux> flux;
    afw::table::Key<meas::base::FluxErrElement> fluxSigma;
    afw::table::Key<afw::table::Flag> fluxFlag;
    afw::table::Key<Scalar> fracDev;
    afw::table::Key<Scalar> objective;
    afw::table::Key<afw::table::Flag> flags[CModelResult::N_FLAGS];
};

} // anonymous

// ------------------- CModelStageData: per-object data we pass around together a lot -----------------------

// Of all the objects related to CModel, CModelStageData is the only that's created per-source (i.e. in the
// apply() methods.  It holds a bunch of stuff related to a single nonlinear fit stage, both things that
// never change (coordinate systems and their transforms) and some things that change a lot (parameters).

namespace {

struct CModelStageData {
    afw::geom::Point2D measSysCenter;       // position of the object in image ("meas") coordinates
    PTR(afw::coord::Coord) position;        // position of the object in ra,dec
    UnitSystem measSys;                     // coordinate systems for the image being measured
    UnitSystem fitSys;                      // coordinate systems for the model parameters
    LocalUnitTransform fitSysToMeasSys;     // coordinate transform from fitSys to measSys
    ndarray::Array<Scalar,1,1> parameters;  // all free parameters (nonlinear + amplitudes)
    ndarray::Array<Scalar,1,1> nonlinear;   // nonlinear parameters (a view into parameters array)
    ndarray::Array<Scalar,1,1> amplitudes;  // linear parameters (a view into parameters array)
    ndarray::Array<Scalar,1,1> fixed;       // fixed parameters (not being fit, still needed to eval model)
    shapelet::MultiShapeletFunction psf;    // multi-shapelet approximation to PSF

    CModelStageData(
        afw::image::Exposure<Pixel> const & exposure,
        Scalar approxFlux, afw::geom::Point2D const & center,
        shapelet::MultiShapeletFunction const & psf_,
        Model const & model
    ) :
        measSysCenter(center), position(exposure.getWcs()->pixelToSky(center)),
        measSys(exposure), fitSys(*position, exposure.getCalib()->getMagnitude(approxFlux)),
        fitSysToMeasSys(*position, fitSys, measSys),
        parameters(ndarray::allocate(model.getNonlinearDim() + model.getAmplitudeDim())),
        nonlinear(parameters[ndarray::view(0, model.getNonlinearDim())]),
        amplitudes(parameters[ndarray::view(model.getNonlinearDim(), parameters.getSize<0>())]),
        fixed(ndarray::allocate(model.getFixedDim())),
        psf(psf_)
    {}

    CModelStageData changeModel(Model const & model) const {
        // If we allowed centroids to vary in some stages and not others, this would resize the parameter
        // arrays and update them accordingly.  For now we just assert that dimensions haven't changed
        // and do a deep-copy.
        // In theory, we should also assert that the ellipse parametrizations haven't changed, but that
        // assert would be too much work to be worthwhile, since at present all Models use the same
        // ellipse parametrization
        assert(model.getNonlinearDim() == nonlinear.getSize<0>());
        assert(model.getAmplitudeDim() == amplitudes.getSize<0>());
        assert(model.getFixedDim() == fixed.getSize<0>());
        CModelStageData r(*this);
        r.parameters = ndarray::copy(parameters);
        r.nonlinear = r.parameters[ndarray::view(0, model.getNonlinearDim())];
        r.amplitudes = r.parameters[ndarray::view(model.getNonlinearDim(), parameters.getSize<0>())];
        // don't need to deep-copy fixed parameters because they're, well, fixed
        return r;
    }

};

} // anonymous

// ------------------- Private Implementation objects -------------------------------------------------------

namespace {

// utility function to create a model matrix: just allocates space for the matrix and calls the likelihood
// object to do the work.
ndarray::Array<Pixel,2,-1> makeModelMatrix(
    Likelihood const & likelihood,
    ndarray::Array<Scalar const,1,1> const & nonlinear
) {
    ndarray::Array<Pixel,2,2> modelMatrixT
        = ndarray::allocate(likelihood.getAmplitudeDim(), likelihood.getDataDim());
    ndarray::Array<Pixel,2,-1> modelMatrix = modelMatrixT.transpose();
    likelihood.computeModelMatrix(modelMatrix, nonlinear);
    return modelMatrix;
}

// Implementation object for a single nonlinear stage (one of "initial", "exp", "dev")
// Note that this doesn't hold its own CModelStageControl; that's held by the CModelControl
// in the main CModelAlgorithm class (for historical and compatibility-with-HSC-fork reasons),
// and hence passed to every method here that needs it.
class CModelStageImpl {
public:
    shapelet::RadialProfile const * profile; // what profile we're trying to fit (ref to singleton)
    PTR(Model) model;                        // defition of parameters, and how to map to Gaussians
    PTR(Prior) prior;                        // Bayesian prior on parameters
    mutable Model::EllipseVector ellipses;   // workspace for asking Model to turn parameters into ellipses
    PTR(afw::table::BaseTable) historyTable;       // optimizer trace Table object
    PTR(OptimizerHistoryRecorder) historyRecorder; // optimizer trace keys/handler

    explicit CModelStageImpl(CModelStageControl const & ctrl) :
        profile(&ctrl.getProfile()),
        model(ctrl.getModel()),
        prior(ctrl.getPrior()),
        ellipses(model->makeEllipseVector())
    {
        if (ctrl.doRecordHistory) {
            afw::table::Schema historySchema;
            historyRecorder.reset(new OptimizerHistoryRecorder(historySchema, model, true));
            historyTable = afw::table::BaseTable::make(historySchema);
        }
    }

    // Create a blank result object, and just fill in the stuff that never changes.
    CModelStageResult makeResult() const {
        CModelStageResult result;
        result.model = model;
        result.prior = prior;
        return result;
    }

    // Use a CModelStageData containing the results of a fit to fill in the higher-level outputs
    // that are part of a CModelStageResult
    void fillResult(
        CModelStageResult & result,
        CModelStageData const & data,
        Scalar amplitudeVariance
    ) const {
        // these are shallow assignments
        result.nonlinear = data.nonlinear;
        result.amplitudes = data.amplitudes;
        result.fixed = data.fixed;
        // flux is just the amplitude converted from fitSys to measSys
        result.flux = data.amplitudes[0] * data.fitSysToMeasSys.flux;
        result.fluxSigma = std::sqrt(amplitudeVariance) * data.fitSysToMeasSys.flux;
        // to compute the ellipse, we need to first read the nonlinear parameters into the workspace
        // ellipse vector, then transform from fitSys to measSys.
        model->writeEllipses(data.nonlinear.begin(), data.fixed.begin(), ellipses.begin());
        result.ellipse = ellipses.front().getCore().transform(data.fitSysToMeasSys.geometric.getLinear());
    }

    // Do the full nonlinear fit for this stage
    void fit(
        CModelStageControl const & ctrl, CModelStageResult & result, CModelStageData const & data,
        afw::image::Exposure<Pixel> const & exposure, afw::detection::Footprint const & footprint
    ) const {
        long long startTime = 0;
        if (ctrl.doRecordTime) {
            startTime = daf::base::DateTime::now().nsecs();
        }
        PTR(UnitTransformedLikelihood) likelihood = boost::make_shared<UnitTransformedLikelihood>(
            model, data.fixed, data.fitSys, *data.position,
            exposure, footprint, data.psf, ctrl.likelihood
        );
        PTR(OptimizerObjective) objective = OptimizerObjective::makeFromLikelihood(likelihood, prior);
        result.objfunc = objective;
        Optimizer optimizer(objective, data.parameters, ctrl.optimizer);
        try {
            if (ctrl.doRecordHistory) {
                result.history = afw::table::BaseCatalog(historyTable);
                optimizer.run(*historyRecorder, result.history);
            } else {
                optimizer.run();
            }
        } catch (std::overflow_error &) {
            result.setFlag(CModelStageResult::NUMERIC_ERROR, true);
        } catch (std::underflow_error &) {
            result.setFlag(CModelStageResult::NUMERIC_ERROR, true);
        } catch (pex::exceptions::UnderflowError &) {
            result.setFlag(CModelStageResult::NUMERIC_ERROR, true);
        } catch (pex::exceptions::OverflowError &) {
            result.setFlag(CModelStageResult::NUMERIC_ERROR, true);
        }

        // Use the optimizer state to set flags.  There's more information in the state than we
        // report in the result, but it's only useful for debugging, and for that the user should
        // look at the history by running outside of plugin mode.
        int state = optimizer.getState();
        if (state & Optimizer::FAILED) {
            result.setFlag(CModelStageResult::FAILED, true);
            if (state & Optimizer::FAILED_MAX_ITERATIONS) {
                result.setFlag(CModelStageResult::MAX_ITERATIONS, true);
            }
        } else {
            result.setFlag(CModelStageResult::FAILED, false);
            if (state & Optimizer::CONVERGED_TR_SMALL) {
                result.setFlag(CModelStageResult::TR_SMALL, true);
            }
        }

        result.objective = optimizer.getObjectiveValue();

        // Set the output parameter vectors.  We deep-assign to the data object to split nonlinear and
        // amplitudes, then shallow-assign these to the result object.
        data.parameters.deep() = optimizer.getParameters(); // sets nonlinear and amplitudes - they are views

        // This amplitudeVariance is computed holding all the nonlinear parameters fixed, which is likely
        // what we'd want for colors, but underestimates the actual uncertainty on the total flux.
        int amplitudeOffset = model->getNonlinearDim();
        Scalar amplitudeVariance = 1.0 / optimizer.getHessian()[amplitudeOffset][amplitudeOffset];

        // Set parameter vectors, flux values, ellipse on result.
        fillResult(result, data, amplitudeVariance);

        if (ctrl.doRecordTime) {
            result.time = (daf::base::DateTime::now().nsecs() - startTime) * 1E9;
        }
    }

    // Do a linear-only fit for this stage (used only in forced mode)
    void fitLinear(
        CModelStageControl const & ctrl, CModelStageResult & result, CModelStageData const & data,
        afw::image::Exposure<Pixel> const & exposure, afw::detection::Footprint const & footprint
    ) const {
        UnitTransformedLikelihood likelihood(
            model, data.fixed, data.fitSys, *data.position,
            exposure, footprint, data.psf, ctrl.likelihood
        );
        ndarray::Array<Pixel,2,-1> modelMatrix = makeModelMatrix(likelihood, data.nonlinear);
        afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(
            modelMatrix,
            likelihood.getData()
        );
        data.amplitudes.deep() = lstsq.getSolution();
        result.objective
            = 0.5*(
                likelihood.getData().asEigen().cast<Scalar>()
                - modelMatrix.asEigen().cast<Scalar>() * lstsq.getSolution().asEigen()
            ).squaredNorm();
        fillResult(result, data, lstsq.getCovariance()[0][0]);
        result.setFlag(CModelStageResult::FAILED, false);
    }

};

} // anonymous

// Master implementation object for CModel.
// Note that this doesn't hold its own CModelControl; that's held by the CModelAlgorithm class
// (for historical and compatibility-with-HSC-fork reasons), and hence passed to every method
// here that needs it.
class CModelAlgorithm::Impl {
public:

    CModelStageImpl initial;  // Implementation object for initial nonlinear fitting stage
    CModelStageImpl exp;      // Implementation object for exponential nonlinear fitting stage
    CModelStageImpl dev;      // Implementation object for de Vaucouleur nonlinear fitting stage
    PTR(Model) model;         // Model object used in final two-component linear fit
    PTR(CModelKeys) keys;     // Key object used to map Result objects to SourceRecord outputs
                              // and extract shapelet PSF approximation.  May be null, depending
                              // on the CModelAlgorithm ctor called
    PTR(CModelKeys) refKeys;  // Key object used to retreive reference ellipses in forced mode
    afw::image::MaskPixel badPixelMask;      // Bitwise OR of mask values we want to leave out of fit region
    std::set<boost::int64_t> diagnosticIds;  // IDs of objects for which we should save diagnostic outputs

    explicit Impl(CModelControl const & ctrl) :
        initial(ctrl.initial), exp(ctrl.exp), dev(ctrl.dev),
        badPixelMask(0x0)
    {
        // turn bad mask plane strings into a bitmask
        for (
            std::vector<std::string>::const_iterator iter = ctrl.region.badMaskPlanes.begin(),
                end = ctrl.region.badMaskPlanes.end();
            iter != end;
            ++iter
        ) {
            badPixelMask |= afw::image::Mask<>::getPlaneBitMask(*iter);
        }

        // construct linear combination model
        ModelVector components(2);
        components[0] = exp.model;
        components[1] = dev.model;
        Model::NameVector prefixes(2);
        prefixes[0] = "exp";
        prefixes[1] = "dev";
        model = boost::make_shared<MultiModel>(components, prefixes);
        // create set of diagnostic IDs for fast lookup
        if (ctrl.diagnostics.enabled) {
            diagnosticIds.insert(ctrl.diagnostics.ids.begin(), ctrl.diagnostics.ids.end());
        }
    }

    // Create a blank result object, filling in only the things that don't change
    CModelResult makeResult() const {
        CModelResult result;
        result.initial = initial.makeResult();
        result.exp = exp.makeResult();
        result.dev = dev.makeResult();
        return result;
    }

    // Do the final two-component linear fit.
    void fitLinear(
        CModelControl const & ctrl, CModelResult & result,
        CModelStageData const & expData, CModelStageData const & devData,
        afw::image::Exposure<Pixel> const & exposure, afw::detection::Footprint const & footprint
    ) const {
        // concatenate exp and dev parameter arrays to make parameter arrays for combined model
        ndarray::Array<Scalar,1,1> nonlinear = ndarray::allocate(model->getNonlinearDim());
        nonlinear[ndarray::view(0, exp.model->getNonlinearDim())] = expData.nonlinear;
        nonlinear[ndarray::view(exp.model->getNonlinearDim(), model->getNonlinearDim())] = devData.nonlinear;
        ndarray::Array<Scalar,1,1> fixed = ndarray::allocate(model->getFixedDim());
        fixed[ndarray::view(0, exp.model->getFixedDim())] = expData.fixed;
        fixed[ndarray::view(exp.model->getFixedDim(), model->getFixedDim())] = devData.fixed;

        UnitTransformedLikelihood likelihood(
            model, fixed, expData.fitSys, *expData.position,
            exposure, footprint, expData.psf, ctrl.likelihood
        );
        ndarray::Array<Pixel,2,-1> modelMatrix = makeModelMatrix(likelihood, nonlinear);
        Vector gradient = -(modelMatrix.asEigen().adjoint() * likelihood.getData().asEigen()).cast<Scalar>();
        Matrix hessian = Matrix::Zero(likelihood.getAmplitudeDim(), likelihood.getAmplitudeDim());
        hessian.selfadjointView<Eigen::Lower>().rankUpdate(modelMatrix.asEigen().adjoint().cast<Scalar>());
        Scalar q0 = 0.5*likelihood.getData().asEigen().squaredNorm();

        // Use truncated Gaussian to compute the maximum-likelihood amplitudes with the constraint
        // that all amplitude must be >= 0
        TruncatedGaussian tg = TruncatedGaussian::fromSeriesParameters(q0, gradient, hessian);
        Vector amplitudes = tg.maximize();
        result.flux = expData.fitSysToMeasSys.flux * amplitudes.sum();

        // To compute the uncertainty on the cmodel flux, we start by transforming the amplitude parameters
        // by the following orthogonal matrix:
        Matrix p(2,2);
        p <<
            M_SQRT1_2, -M_SQRT1_2,
            M_SQRT1_2, M_SQRT1_2;
        Vector p_mu = p * amplitudes;
        Matrix p_hessian_pt = p * hessian * p.adjoint();
        // After this transformation, \sqrt(2)*p_mu[1] is the total flux, and \sqrt(2)*p_mu[0] is the
        // difference between the fluxes of the two components.
        // We define the flux error as the variance on the total flux with the difference between the
        // fluxes held fixed.  This is artificial, and an underestimate of the true uncertainty, as we
        // ought to be marginalizing over the difference between the fluxes - but that integral would
        // often diverge if we don't take into account the constraints, and if we do, the probability
        // distribution we get is complicated (Gaussian times a difference of error functions) and
        // asymmetric, so we'll leave that as a potential project for the future.
        result.fluxSigma = expData.fitSysToMeasSys.flux / std::sqrt(p_hessian_pt(1,1));
        result.setFlag(CModelResult::FAILED, false);

        result.fracDev = amplitudes[1] / amplitudes.sum();
        result.objective = tg.evaluateLog()(amplitudes);
    }

    // Guess parameters for the initial fit stage from image moments
    void guessParametersFromMoments(
        CModelControl const & ctrl, CModelStageData & data,
        afw::geom::ellipses::Quadrupole const & moments
    ) const {

        // Deconvolve the moments ellipse, with a floor to keep the result from
        // having moments <= 0
        afw::geom::ellipses::Ellipse psfEllipse = data.psf.evaluate().computeMoments();
        afw::geom::ellipses::Quadrupole psfMoments(psfEllipse.getCore());
        Scalar mir2 = ctrl.minInitialRadius * ctrl.minInitialRadius;
        Scalar ixx = std::max(moments.getIxx() - psfMoments.getIxx(), mir2);
        Scalar iyy = std::max(moments.getIyy() - psfMoments.getIyy(), mir2);
        Scalar ixy = moments.getIxy() - psfMoments.getIxy();
        if (ixx*iyy < ixy*ixy) {
            ixy = 0.0;
        }
        afw::geom::ellipses::Quadrupole deconvolvedMoments(
            ixx, iyy, ixy,
            true // throw if ellipse is invalid
        );
        afw::geom::ellipses::Ellipse deconvolvedEllipse(
            deconvolvedMoments,
            afw::geom::Point2D(data.measSysCenter - psfEllipse.getCenter())
        );

        // Convert ellipse from moments to half-light using the ratio for this profile
        deconvolvedEllipse.getCore().scale(1.0 / initial.profile->getMomentsRadiusFactor());

        // Transform the deconvolved ellipse from MeasSys to FitSys
        deconvolvedEllipse.transform(data.fitSysToMeasSys.geometric.invert()).inPlace();

        // Convert to the ellipse parametrization used by the Model (assigning to an ellipse converts
        // between parametrizations)
        assert(initial.ellipses.size() == 1u); // should be true of all Models that come from RadialProfiles
        initial.ellipses.front() = deconvolvedEllipse;

        // Read the ellipse into the nonlinear and fixed parameters.
        initial.model->readEllipses(initial.ellipses.begin(), data.nonlinear.begin(), data.fixed.begin());

        // Set the initial amplitude (a.k.a. flux) to 1: recall that in FitSys, this is approximately correct
        assert(data.amplitudes.getSize<0>() == 1); // should be true of all Models from RadialProfiles
        data.amplitudes[0] = 1.0;

        // Ensure the initial parameters are compatible with the prior
        if (initial.prior && initial.prior->evaluate(data.nonlinear, data.amplitudes) == 0.0) {
            initial.ellipses.front().setCore(afw::geom::ellipses::Quadrupole(mir2, mir2, 0.0));
            initial.model->readEllipses(initial.ellipses.begin(), data.nonlinear.begin(), data.fixed.begin());
            if (initial.prior->evaluate(data.nonlinear, data.amplitudes) == 0.0) {
                throw LSST_EXCEPT(
                    meas::base::FatalAlgorithmError,
                    "minInitialRadius is incompatible with prior"
                );
            }
        }
    }

    // Save diagnostic  outputs for a particular source
    template <typename T>
    void writeDiagnostics(
        CModelControl const & ctrl,
        boost::int64_t id,
        CModelResult const & result,
        afw::image::Exposure<T> const & exposure
    ) const {
        if (!result.initialFitRegion) {
            return; // cannot write diagnostics if we didn't at least get this far.
        }
        std::string path = (boost::format("%s/%d.fits") % ctrl.diagnostics.root % id).str();
        afw::fits::Fits fits(path, "w", afw::fits::Fits::AUTO_CLOSE | afw::fits::Fits::AUTO_CHECK);
        afw::geom::Box2I bbox = result.initialFitRegion->getBBox();
        if (result.finalFitRegion) {
            bbox.include(result.finalFitRegion->getBBox());
        }
        afw::image::Image<T> subImage(*exposure.getMaskedImage().getImage(), bbox, afw::image::PARENT);
        subImage.writeFits(fits);
        assert(fits.countHdus() == 1);
        if (ctrl.initial.doRecordHistory && result.initial.history.getTable()) {
            result.initial.history.writeFits(fits);
        } else {
            fits.createEmpty();
        }
        assert(fits.countHdus() == 2);
        if (ctrl.exp.doRecordHistory && result.exp.history.getTable()) {
            result.exp.history.writeFits(fits);
        } else {
            fits.createEmpty();
        }
        assert(fits.countHdus() == 3);
        if (ctrl.dev.doRecordHistory && result.dev.history.getTable()) {
            result.dev.history.writeFits(fits);
        } else {
            fits.createEmpty();
        }
        assert(fits.countHdus() == 4);
    }

};

// ------------------- CModelAlgorithm itself ---------------------------------------------------------------

CModelAlgorithm::CModelAlgorithm(
    std::string const & name,
    Control const & ctrl,
    afw::table::Schema & schema
) : _ctrl(ctrl), _impl(new Impl(ctrl))
{
    _impl->keys = boost::make_shared<CModelKeys>(
        *_impl->initial.model, *_impl->exp.model, *_impl->dev.model,
        boost::ref(schema), name, false, ctrl
    );
}

CModelAlgorithm::CModelAlgorithm(
    std::string const & name,
    Control const & ctrl,
    afw::table::SchemaMapper & schemaMapper
) : _ctrl(ctrl), _impl(new Impl(ctrl))
{
    _impl->keys = boost::make_shared<CModelKeys>(
        *_impl->initial.model, *_impl->exp.model, *_impl->dev.model,
        boost::ref(schemaMapper.editOutputSchema()), name, true, ctrl
    );
    _impl->refKeys = boost::make_shared<CModelKeys>(
        *_impl->initial.model, *_impl->exp.model, *_impl->dev.model,
        schemaMapper.getInputSchema(), name
    );
}

CModelAlgorithm::CModelAlgorithm(Control const & ctrl) :
    _ctrl(ctrl), _impl(new Impl(ctrl))
{}

PTR(afw::detection::Footprint) CModelAlgorithm::determineInitialFitRegion(
    afw::image::Mask<> const & mask,
    afw::detection::Footprint const & footprint,
    afw::geom::Box2I const & psfBBox
) const {
    PTR(afw::detection::Footprint) region;
    if (footprint.getArea() > getControl().region.maxArea) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            "Maximum area exceeded by original footprint",
            CModelResult::MAX_AREA
        );
    }
    region = afw::detection::growFootprint(
        footprint,
        getControl().region.nGrowFootprint,
        true
    );
    if (region->getArea() > getControl().region.maxArea) {
        throw LSST_EXCEPT(
             meas::base::MeasurementError,
             "Maximum area exceeded by grown footprint",
             CModelResult::MAX_AREA
        );
    }
    if (getControl().region.includePsfBBox && !region->getBBox().contains(psfBBox)) {
        region = _mergeFootprints(*region, afw::detection::Footprint(psfBBox));
    }
    int originalArea = region->getArea();
    region->clipTo(mask.getBBox(afw::image::PARENT));
    region->intersectMask(mask, _impl->badPixelMask);
    if (originalArea - region->getArea() > originalArea*getControl().region.maxBadPixelFraction) {
        throw LSST_EXCEPT(
             meas::base::MeasurementError,
             "Fraction of bad pixels in fit region exceeds threshold",
             CModelResult::MAX_BAD_PIXEL_FRACTION
        );
    }
    if (region->getArea() > getControl().region.maxArea) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            "Maximum area exceeded by initial fit region",
            CModelResult::MAX_AREA
        );
    }
    return region;
}

PTR(afw::detection::Footprint) CModelAlgorithm::determineFinalFitRegion(
    afw::image::Mask<> const & mask,
    afw::detection::Footprint const & footprint,
    afw::geom::Box2I const & psfBBox,
    afw::geom::ellipses::Ellipse const & ellipse
) const {
    PTR(afw::detection::Footprint) region;
    if (footprint.getArea() > getControl().region.maxArea) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            "Maximum area exceeded by original footprint",
            CModelResult::MAX_AREA
        );
    }
    region = afw::detection::growFootprint(
        footprint,
        getControl().region.nGrowFootprint,
        true
    );
    afw::geom::ellipses::Ellipse fullEllipse(ellipse);
    fullEllipse.getCore().scale(getControl().region.nInitialRadii);
    if (fullEllipse.getCore().getArea() > getControl().region.maxArea) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            "Maximum area exceeded by ellipse component of region",
            CModelResult::MAX_AREA
        );
    }
    afw::detection::Footprint ellipseFootprint(fullEllipse);
    if (ellipseFootprint.getArea() > 0) {
        region = _mergeFootprints(*region, ellipseFootprint);
    }
    if (getControl().region.includePsfBBox && !region->getBBox().contains(psfBBox)) {
        region = _mergeFootprints(*region, afw::detection::Footprint(psfBBox));
    }
    double originalArea = region->getArea();
    region->clipTo(mask.getBBox(afw::image::PARENT));
    region->intersectMask(mask, _impl->badPixelMask);
    if ((1.0 - region->getArea() / originalArea) > getControl().region.maxBadPixelFraction) {
        throw LSST_EXCEPT(
             meas::base::MeasurementError,
             "Fraction of bad pixels in fit region exceeds threshold",
             CModelResult::MAX_BAD_PIXEL_FRACTION
        );
    }
    if (region->getArea() > getControl().region.maxArea) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            "Maximum area exceeded by final fit region",
            CModelResult::MAX_AREA
        );
    }
    return region;
}

CModelAlgorithm::Result CModelAlgorithm::apply(
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    afw::geom::Point2D const & center,
    afw::geom::ellipses::Quadrupole const & moments,
    Scalar approxFlux
) const {
    Result result = _impl->makeResult();
    _applyImpl(result, exposure, footprint, psf, center, moments, approxFlux);
    return result;
}

void CModelAlgorithm::_applyImpl(
    Result & result,
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    afw::geom::Point2D const & center,
    afw::geom::ellipses::Quadrupole const & moments,
    Scalar approxFlux
) const {

    afw::geom::Box2I psfBBox = exposure.getPsf()->computeImage(center)->getBBox(afw::image::PARENT);

    // Grow the footprint, clip bad pixels and the exposure bbox
    result.initialFitRegion = determineInitialFitRegion(
        *exposure.getMaskedImage().getMask(),
        footprint,
        psfBBox
    );

    // Negative approxFlux means we should come up with an estimate ourselves.
    // This is only used to avoid scaling problems in the optimizer, so it doesn't have to be very good.
    if (approxFlux < 0.0) {
        approxFlux = computeFluxInFootprint(*exposure.getMaskedImage().getImage(), footprint);
    }

    // Set up coordinate systems and empty parameter vectors
    CModelStageData initialData(exposure, approxFlux, center, psf, *_impl->initial.model);

    // Initialize the parameter vectors by doing deconvolving the moments
    _impl->guessParametersFromMoments(getControl(), initialData, moments);

    // Do the initial fit
    // TODO: use only 0th-order terms in psf
    _impl->initial.fit(getControl().initial, result.initial, initialData, exposure, *result.initialFitRegion);

    if (result.initial.getFlag(CModelStageResult::FAILED)) return;

    // Include a multiple of the initial-fit ellipse in the footprint, re-do clipping
    result.initial.model->writeEllipses(initialData.nonlinear.begin(), initialData.fixed.begin(),
                                        _impl->initial.ellipses.begin());
    _impl->initial.ellipses.front().transform(initialData.fitSysToMeasSys.geometric).inPlace();

    // Revisit the pixel region to use in the fit, taking into account the initial ellipse
    result.finalFitRegion = determineFinalFitRegion(
        *exposure.getMaskedImage().getMask(),
        footprint,
        psfBBox,
        _impl->initial.ellipses.front()
    );

    // Do the exponential fit
    CModelStageData expData = initialData.changeModel(*_impl->exp.model);
    _impl->exp.fit(getControl().exp, result.exp, expData, exposure, *result.finalFitRegion);

    // Do the de Vaucouleur fit
    CModelStageData devData = initialData.changeModel(*_impl->dev.model);
    _impl->dev.fit(getControl().dev, result.dev, devData, exposure, *result.finalFitRegion);

    if (result.exp.getFlag(CModelStageResult::FAILED) ||result.dev.getFlag(CModelStageResult::FAILED))
        return;

    // Do the linear combination fit
    _impl->fitLinear(getControl(), result, expData, devData, exposure, *result.finalFitRegion);
}

CModelAlgorithm::Result CModelAlgorithm::applyForced(
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    afw::geom::Point2D const & center,
    CModelResult const & reference,
    Scalar approxFlux
) const {
    Result result = _impl->makeResult();
    _applyForcedImpl(result, exposure, footprint, psf, center, reference, approxFlux);
    return result;
}

void CModelAlgorithm::writeResultToRecord(
    Result const & result,
    afw::table::BaseRecord & record
) const {
    if (!_impl->keys) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "Algorithm was not initialized with a schema; cannot copy to record"
        );
    }
    _impl->keys->copyResultToRecord(result, record);
}

void CModelAlgorithm::fail(
    afw::table::SourceRecord & record,
    meas::base::MeasurementError * error
) const {
    if (!_impl->keys) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "Algorithm was not initialized with a schema; cannot handle failures"
        );
    }
    record.set(_impl->keys->flags[CModelResult::FAILED], true);
    if (error) {
        record.set(_impl->keys->flags[error->getFlagBit()], true);
    }
}

void CModelAlgorithm::_applyForcedImpl(
    Result & result,
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    afw::geom::Point2D const & center,
    CModelResult const & reference,
    Scalar approxFlux
) const {
    afw::geom::Box2I psfBBox = exposure.getPsf()->computeImage(center)->getBBox(afw::image::PARENT);

    // Negative approxFlux means we should come up with an estimate ourselves.
    // This is only used to avoid scaling problems in the optimizer, so it doesn't have to be very good.
    if (approxFlux < 0.0) {
        approxFlux = computeFluxInFootprint(*exposure.getMaskedImage().getImage(), footprint);
    }

    // Set up coordinate systems and empty parameter vectors
    CModelStageData initialData(exposure, approxFlux, center, psf, *_impl->initial.model);

    // Initialize the parameter vectors from the reference values.  Because these are
    // in fitSys units, we don't need to transform them, as fitSys (or at least its
    // Wcs) should be the same in both forced mode and non-forced mode.
    initialData.nonlinear.deep() = reference.initial.nonlinear;
    initialData.fixed.deep() = reference.initial.fixed;
    // Read those parameters into the ellipses.
    _impl->initial.model->writeEllipses(initialData.nonlinear.begin(), initialData.fixed.begin(),
                                        _impl->initial.ellipses.begin());
    // Transform the ellipses to the exposure coordinate system
    _impl->initial.ellipses.front().transform(initialData.fitSysToMeasSys.geometric).inPlace();

    // Grow the footprint and include the initial ellipse, clip bad pixels and the exposure bbox;
    // in forced mode we can just use the final fit region immediately since we won't be changing
    // the initial fit ellipse.
    result.finalFitRegion = determineFinalFitRegion(
        *exposure.getMaskedImage().getMask(),
        footprint,
        psfBBox,
        _impl->initial.ellipses.front()
    );

    // Do the initial fit (amplitudes only)
    if (!reference.initial.getFlag(CModelStageResult::FAILED)) {
        _impl->initial.fitLinear(getControl().initial, result.initial, initialData,
                                 exposure, *result.finalFitRegion);
    }

    // Do the exponential fit (amplitudes only)
    CModelStageData expData = initialData.changeModel(*_impl->exp.model);
    if (!reference.exp.getFlag(CModelStageResult::FAILED)) {
        expData.nonlinear.deep() = reference.exp.nonlinear;
        expData.fixed.deep() = reference.exp.fixed;
        _impl->exp.fitLinear(getControl().exp, result.exp, expData, exposure, *result.finalFitRegion);
    }

    // Do the de Vaucouleur fit (amplitudes only)
    CModelStageData devData = initialData.changeModel(*_impl->dev.model);
    if (!reference.dev.getFlag(CModelStageResult::FAILED)) {
        devData.nonlinear.deep() = reference.dev.nonlinear;
        devData.fixed.deep() = reference.dev.fixed;
        _impl->dev.fitLinear(getControl().dev, result.dev, devData, exposure, *result.finalFitRegion);
    }

    if (result.exp.getFlag(CModelStageResult::FAILED) ||result.dev.getFlag(CModelStageResult::FAILED))
        return;

    // Do the linear combination fit
    _impl->fitLinear(getControl(), result, expData, devData, exposure, *result.finalFitRegion);
}

template <typename PixelT>
shapelet::MultiShapeletFunction CModelAlgorithm::_processInputs(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure
) const {
    // Set all failure flags so that's the result if we throw.
    source.set(_impl->keys->flags[Result::FAILED], true);
    source.set(_impl->keys->initial.flags[CModelStageResult::FAILED], true);
    source.set(_impl->keys->exp.flags[CModelStageResult::FAILED], true);
    source.set(_impl->keys->dev.flags[CModelStageResult::FAILED], true);
    if (!_impl->keys) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "Algorithm was not initialized with a schema; cannot run in plugin mode"
        );
    }
    if (!exposure.getWcs()) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "Exposure has no Wcs"
        );
    }
    if (!exposure.getCalib() || exposure.getCalib()->getFluxMag0().first == 0.0) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "Exposure has no valid Calib"
        );
    }
    if (!exposure.getPsf()) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "Exposure has no Psf"
        );
    }
    return source.get(_impl->keys->psf);
}

void CModelAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<Pixel> const & exposure
) const {
    Result result = _impl->makeResult();
    // Read the shapelet approximation to the PSF, load/verify other inputs from the SourceRecord
    shapelet::MultiShapeletFunction psf = _processInputs(measRecord, exposure);
    if (!measRecord.getTable()->getShapeKey().isValid() ||
        (measRecord.getTable()->getShapeFlagKey().isValid() && measRecord.getShapeFlag())) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            "Shape slot algorithm failed or was not run",
            Result::NO_SHAPE
        );
    }
    // If PsfFlux has been run, use that for approx flux; otherwise we'll compute it ourselves.
    Scalar approxFlux = -1.0;
    if (measRecord.getTable()->getPsfFluxKey().isValid() && !measRecord.getPsfFluxFlag()) {
        approxFlux = measRecord.getPsfFlux();
    }
    try {
        _applyImpl(result, exposure, *measRecord.getFootprint(), psf, measRecord.getCentroid(),
                   measRecord.getShape(), approxFlux);
    } catch (...) {
        _impl->keys->copyResultToRecord(result, measRecord);
        if (_impl->diagnosticIds.find(measRecord.getId()) != _impl->diagnosticIds.end()) {
            _impl->writeDiagnostics(getControl(), measRecord.getId(), result, exposure);
        }
        throw;
    }
    _impl->keys->copyResultToRecord(result, measRecord);
    if (_impl->diagnosticIds.find(measRecord.getId()) != _impl->diagnosticIds.end()) {
        _impl->writeDiagnostics(getControl(), measRecord.getId(), result, exposure);
    }
}

void CModelAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<Pixel> const & exposure,
    afw::table::SourceRecord const & refRecord
) const {
    Result result = _impl->makeResult();
    assert(measRecord.getFootprint()->getArea());
    // Read the shapelet approximation to the PSF, load/verify other inputs from the SourceRecord
    shapelet::MultiShapeletFunction psf = _processInputs(measRecord, exposure);
    // If PsfFlux has been run, use that for approx flux; otherwise we'll compute it ourselves.
    Scalar approxFlux = -1.0;
    if (measRecord.getTable()->getPsfFluxKey().isValid() && !measRecord.getPsfFluxFlag()) {
        approxFlux = measRecord.getPsfFlux();
    }
    try {
        Result refResult = _impl->refKeys->copyRecordToResult(refRecord);
        _applyForcedImpl(result, exposure, *measRecord.getFootprint(), psf, measRecord.getCentroid(),
                         refResult, approxFlux);
    } catch (...) {
        _impl->keys->copyResultToRecord(result, measRecord);
        throw;
    }
    _impl->keys->copyResultToRecord(result, measRecord);
}

}}} // namespace lsst::meas::multifit
