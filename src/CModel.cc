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
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/meas/extensions/multiShapelet/FitPsf.h"
#include "lsst/meas/multifit/TruncatedGaussian.h"
#include "lsst/meas/multifit/CModel.h"

namespace lsst { namespace meas { namespace multifit {


//-------------------- Utility code -------------------------------------------------------------------------

namespace {

PTR(afw::detection::Footprint) mergeFootprints(
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
            pex::exceptions::RuntimeErrorException,
            "Footprints to be merged do not overlap"
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
    char const * pkgDir = std::getenv("MEAS_MULTIFIT_DIR");
    if (!pkgDir) {
        throw LSST_EXCEPT(
            pex::exceptions::IoErrorException,
            "MEAS_MULTIFIT_DIR environment variable not defined; cannot find persisted Priors"
        );
    }
    boost::filesystem::path priorPath
        = boost::filesystem::path(pkgDir)
        / boost::filesystem::path("data")
        / boost::filesystem::path(priorName + ".fits");
    PTR(Mixture) mixture = Mixture::readFits(priorPath.string());
    return boost::make_shared<MixturePrior>(mixture, "single-ellipse");
}

PTR(algorithms::AlgorithmControl) CModelControl::_clone() const {
    return boost::make_shared<CModelControl>(*this);
}

PTR(algorithms::Algorithm) CModelControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    return boost::make_shared<CModelAlgorithm>(*this, boost::ref(schema), others, isForced);
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

struct CModelStageKeys {

    // this constructor is used to allocate output fields in both forced and non-forced mode
    CModelStageKeys(
        Model const & model,
        afw::table::Schema & schema,
        std::string const & prefix,
        std::string const & stage,
        bool isForced
    ) :
        flux(afw::table::addFluxFields(schema, prefix + ".flux", "flux from the " + stage + " fit")),
        fluxCorrection(prefix, schema)
    {
        if (!isForced) {
            ellipse = schema.addField<afw::table::Moments<Scalar> >(
                prefix + ".ellipse", "effective radius ellipse from the " + stage + " fit"
            );
            objective = schema.addField<Scalar>(
                prefix + ".objective", "-ln(likelihood*prior) at best-fit point for the " + stage + " fit"
            );
            nonlinear = schema.addField<afw::table::Array<Scalar> >(
                prefix + ".nonlinear", "nonlinear parameters for the " + stage + " fit",
                model.getNonlinearDim()
            );
            fixed = schema.addField<afw::table::Array<Scalar> >(
                prefix + ".fixed", "fixed parameters for the " + stage + " fit",
                model.getFixedDim()
            );
            flags[CModelStageResult::TR_SMALL] = schema.addField<afw::table::Flag>(
                prefix + ".flags.trSmall",
                "the optimizer converged because the trust radius became too small; this is a less-secure "
                "result than when the gradient is below the threshold, but usually not a problem"
            );
            flags[CModelStageResult::MAX_ITERATIONS] = schema.addField<afw::table::Flag>(
                prefix + ".flags.maxIter",
                "the optimizer hit the maximum number of iterations and did not converge"
            );
        }
        flags[CModelStageResult::FAILED] = flux.flag; // these flags refer to the same underlying field
    }

    // this constructor is used to get needed keys from the reference schema in forced mode
    CModelStageKeys(
        Model const & model,
        afw::table::Schema const & schema,
        std::string const & prefix
    ) :
        flux(schema[prefix + ".flux"], schema[prefix + ".flux.err"], schema[prefix + ".flags"]),
        nonlinear(schema[prefix + ".nonlinear"]),
        fixed(schema[prefix + ".fixed"])
    {
        flags[CModelStageResult::FAILED] = flux.flag; // these flags refer to the same underlying field
        LSST_THROW_IF_NE(
            model.getNonlinearDim(), nonlinear.getSize(),
            pex::exceptions::LengthErrorException,
            "Configured model nonlinear dimension (%d) does not match reference schema (%d)"
        );
        LSST_THROW_IF_NE(
            model.getFixedDim(), fixed.getSize(),
            pex::exceptions::LengthErrorException,
            "Configured model fixed dimension (%d) does not match reference schema (%d)"
        );
    }

    void copyResultToRecord(CModelStageResult const & result, afw::table::BaseRecord & record) {
        record.set(flux.meas, result.flux);
        record.set(flux.err, result.fluxSigma);
        record.set(flux.flag, result.getFlag(CModelStageResult::FAILED));
        record.set(fluxCorrection.psfFactor, 1.0); // TODO
        record.set(fluxCorrection.psfFactorFlag, false); // TODO
        record.set(objective, result.objective);
        if (ellipse.isValid()) {
            record.set(ellipse, result.ellipse);
        }
        if (nonlinear.isValid() && !result.nonlinear.isEmpty()) {
            record.set(nonlinear, result.nonlinear);
        }
        if (fixed.isValid() && !result.fixed.isEmpty()) {
            record.set(fixed, result.fixed);
        }
        for (int b = 0; b < CModelStageResult::N_FLAGS; ++b) {
            record.set(flags[b], result.flags[b]);
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

    afw::table::KeyTuple<afw::table::Flux> flux;
    algorithms::ScaledFlux::KeyTuple fluxCorrection;
    afw::table::Key<afw::table::Moments<Scalar> > ellipse;
    afw::table::Key<Scalar> objective;
    afw::table::Key<afw::table::Flag> flags[CModelStageResult::N_FLAGS];
    afw::table::Key<afw::table::Array<Scalar> > nonlinear;
    afw::table::Key<afw::table::Array<Scalar> > fixed;
};

struct CModelKeys {

    // this constructor is used to allocate output fields in both forced and non-forced mode
    CModelKeys(
        Model const & initialModel, Model const & expModel, Model const & devModel,
        afw::table::Schema & schema,
        std::string const & prefix,
        bool isForced
    ) :
        initial(initialModel, schema, prefix + ".initial", "initial", isForced),
        exp(expModel, schema, prefix + ".exp", "exponential", isForced),
        dev(devModel, schema, prefix + ".dev", "de Vaucouleur", isForced),
        flux(afw::table::addFluxFields(schema, prefix + ".flux", "flux from the final cmodel fit")),
        fluxCorrection(prefix, schema),
        fracDev(schema.addField<Scalar>(prefix + ".fracDev", "fraction of flux in de Vaucouleur component")),
        objective(schema.addField<Scalar>(prefix + ".objective", "-ln(likelihood) (chi^2) in cmodel fit"))
    {
        flags[CModelResult::FAILED] = flux.flag; // these keys refer to the same underlying field
        flags[CModelResult::MAX_AREA] = schema.addField<afw::table::Flag>(
            prefix + ".flags.maxArea",
            "number of pixels in fit region exceeded the region.maxArea value (usually due to bad moments)"
        );
        flags[CModelResult::MAX_BAD_PIXEL_FRACTION] = schema.addField<afw::table::Flag>(
            prefix + ".flags.maxBadPixelFraction",
            "the fraction of bad/clipped pixels in the fit region exceeded region.maxBadPixelFraction"
        );
        flags[CModelResult::NO_SHAPE] = schema.addField<afw::table::Flag>(
            prefix + ".flags.noShape",
            "the shape slot needed to initialize the parameters failed or was not defined"
        );
        flags[CModelResult::NO_PSF] = schema.addField<afw::table::Flag>(
            prefix + ".flags.noPsf",
            "the multishapelet fit to the PSF model did not succeed"
        );
        flags[CModelResult::NO_WCS] = schema.addField<afw::table::Flag>(
            prefix + ".flags.noWcs",
            "input exposure has no world coordinate system information"
        );
        flags[CModelResult::NO_CALIB] = schema.addField<afw::table::Flag>(
            prefix + ".flags.noCalib",
            "input exposure has no photometric calibration information"
        );
    }

    // this constructor is used to get needed keys from the reference schema in forced mode
    CModelKeys(
        Model const & initialModel, Model const & expModel, Model const & devModel,
        afw::table::Schema const & schema,
        std::string const & prefix
    ) :
        initial(initialModel, schema, prefix + ".initial"),
        exp(expModel, schema, prefix + ".exp"),
        dev(devModel, schema, prefix + ".dev")
    {}

    void copyResultToRecord(CModelResult const & result, afw::table::BaseRecord & record) {
        initial.copyResultToRecord(result.initial, record);
        exp.copyResultToRecord(result.exp, record);
        dev.copyResultToRecord(result.dev, record);
        record.set(flux.meas, result.flux);
        record.set(flux.err, result.fluxSigma);
        record.set(fluxCorrection.psfFactor, 1.0); // TODO
        record.set(fluxCorrection.psfFactorFlag, false); // TODO
        record.set(fracDev, result.fracDev);
        record.set(objective, result.objective);
        for (int b = 0; b < CModelStageResult::N_FLAGS; ++b) {
            record.set(flags[b], result.flags[b]);
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
    afw::table::KeyTuple<afw::table::Flux> flux;
    algorithms::ScaledFlux::KeyTuple fluxCorrection;
    afw::table::Key<Scalar> fracDev;
    afw::table::Key<Scalar> objective;
    afw::table::Key<afw::table::Flag> flags[CModelResult::N_FLAGS];
};

} // anonymous

// ------------------- CModelStageData: per-object data we pass around together a lot -----------------------

namespace {

struct CModelStageData {
    afw::geom::Point2D measSysCenter;
    PTR(afw::coord::Coord) position;
    UnitSystem measSys;
    UnitSystem fitSys;
    LocalUnitTransform fitSysToMeasSys;
    ndarray::Array<Scalar,1,1> parameters;
    ndarray::Array<Scalar,1,1> nonlinear;
    ndarray::Array<Scalar,1,1> amplitudes;
    ndarray::Array<Scalar,1,1> fixed;
    shapelet::MultiShapeletFunction psf;

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
        // If we allow centroids to vary in some stages and not others, this will resize the parameter
        // arrays and update them accordingly.  For now we just assert that dimensions haven't changed
        // and do a deep-copy.
        // In theory, we should also assert that the ellipse parametrizations haven't changed, but that
        // assert would be too much work to be worthwhile.
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

class CModelStageImpl {
public:
    shapelet::RadialProfile const * profile;
    PTR(Model) model;
    PTR(Prior) prior;
    mutable Model::EllipseVector ellipses;
    PTR(afw::table::BaseTable) historyTable;
    PTR(OptimizerHistoryRecorder) historyRecorder;

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

    CModelStageResult makeResult() const {
        CModelStageResult result;
        result.model = model;
        result.prior = prior;
        return result;
    }

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

    void fit(
        CModelStageControl const & ctrl, CModelStageResult & result, CModelStageData const & data,
        afw::image::Exposure<Pixel> const & exposure, afw::detection::Footprint const & footprint
    ) const {
        PTR(ProjectedLikelihood) likelihood = boost::make_shared<ProjectedLikelihood>(
            model, data.fixed, data.fitSys, *data.position,
            exposure, footprint, data.psf, ctrl.likelihood
        );
        PTR(OptimizerObjective) objective = OptimizerObjective::makeFromLikelihood(likelihood, prior);
        Optimizer optimizer(objective, data.parameters, ctrl.optimizer);
        if (ctrl.doRecordHistory) {
            result.history = afw::table::BaseCatalog(historyTable);
            optimizer.run(*historyRecorder, result.history);
        } else {
            optimizer.run();
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
    }

    void fitLinear(
        CModelStageControl const & ctrl, CModelStageResult & result, CModelStageData const & data,
        afw::image::Exposure<Pixel> const & exposure, afw::detection::Footprint const & footprint
    ) const {
        ProjectedLikelihood likelihood(
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


class CModelAlgorithm::Impl {
public:

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
    }

    CModelStageImpl initial;
    CModelStageImpl exp;
    CModelStageImpl dev;
    PTR(Model) model;
    PTR(CModelKeys) keys;
    PTR(CModelKeys) refKeys;
    PTR(extensions::multiShapelet::FitPsfControl const) fitPsfCtrl;
    afw::image::MaskPixel badPixelMask;

    CModelResult makeResult() const {
        CModelResult result;
        result.initial = initial.makeResult();
        result.exp = exp.makeResult();
        result.dev = dev.makeResult();
        return result;
    }

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

        ProjectedLikelihood likelihood(
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

    void guessParametersFromMoments(
        CModelControl const & ctrl, CModelStageData & data,
        afw::geom::ellipses::Quadrupole const & moments
    ) const {

        // Deconvolve the moments ellipse, with a floor to keep the result from
        // having moments <= 0
        afw::geom::ellipses::Ellipse psfEllipse = data.psf.evaluate().computeMoments();
        afw::geom::ellipses::Quadrupole psfMoments(psfEllipse.getCore());
        afw::geom::ellipses::Quadrupole deconvolvedMoments(
            std::max(moments.getIxx() - psfMoments.getIxx(), ctrl.minInitialRadius),
            std::max(moments.getIyy() - psfMoments.getIyy(), ctrl.minInitialRadius),
            moments.getIxy() - psfMoments.getIxy(),
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
    }

};

// ------------------- CModelAlgorithm itself ---------------------------------------------------------------

CModelAlgorithm::CModelAlgorithm(
    Control const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmMap const & others,
    bool isForced
) : algorithms::Algorithm(ctrl), _impl(new Impl(ctrl))
{
    _impl->keys = boost::make_shared<CModelKeys>(
        *_impl->initial.model, *_impl->exp.model, *_impl->dev.model, boost::ref(schema), ctrl.name, isForced
    );
    // Ideally we'd like to initalize refKeys here too when isForced==true, but we aren't passed the
    // refSchema here, so instead we'll construct that on first use.  This will be fixed in the next
    // version of the measurement framework that's in progress on the LSST side.

    algorithms::AlgorithmMap::const_iterator i = others.find(ctrl.psfName);
    if (i == others.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("FitPsf with name '%s' not found; needed by CModel.") % ctrl.psfName).str()
        );
    }
    _impl->fitPsfCtrl = boost::dynamic_pointer_cast<extensions::multiShapelet::FitPsfControl const>(
        i->second->getControl().clone()
    );
    if (!_impl->fitPsfCtrl) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("Algorithm with name '%s' is not FitPsf.") % ctrl.psfName).str()
        );
    }
}

CModelAlgorithm::CModelAlgorithm(Control const & ctrl) :
    algorithms::Algorithm(ctrl), _impl(new Impl(ctrl))
{}

int CModelAlgorithm::getFluxCount() const { return 4; }

afw::table::KeyTuple<afw::table::Flux> CModelAlgorithm::getFluxKeys(int n) const {
    if (!_impl->keys) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Algorithm was not initialized with a schema; no keys allocated"
        );
    }
    switch (n) {
    case 0:
        return _impl->keys->flux;
    case 1:
        return _impl->keys->initial.flux;
    case 2:
        return _impl->keys->exp.flux;
    case 3:
        return _impl->keys->dev.flux;
    default:
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Invalid index for getFluxKeys()"
        );
    }
}

algorithms::ScaledFlux::KeyTuple CModelAlgorithm::getFluxCorrectionKeys(int n) const {
    if (!_impl->keys) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Algorithm was not initialized with a schema; no keys allocated"
        );
    }
    switch (n) {
    case 0:
        return _impl->keys->fluxCorrection;
    case 1:
        return _impl->keys->initial.fluxCorrection;
    case 2:
        return _impl->keys->exp.fluxCorrection;
    case 3:
        return _impl->keys->dev.fluxCorrection;
    default:
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Invalid index for getFluxCorrectionKeys()"
        );
    }
}

PTR(afw::detection::Footprint) CModelAlgorithm::determineInitialFitRegion(
    afw::image::Mask<> const & mask,
    afw::detection::Footprint const & footprint
) const {
    PTR(afw::detection::Footprint) region = afw::detection::growFootprint(
        footprint,
        getControl().region.nGrowFootprint,
        true
    );
    double originalArea = region->getArea();
    region->clipTo(mask.getBBox(afw::image::PARENT));
    region->intersectMask(mask, _impl->badPixelMask);
    if ((1.0 - region->getArea() / originalArea) > getControl().region.maxBadPixelFraction) {
        region.reset();
    }
    return region;
}

PTR(afw::detection::Footprint) CModelAlgorithm::determineFinalFitRegion(
    afw::image::Mask<> const & mask,
    afw::detection::Footprint const & footprint,
    afw::geom::ellipses::Quadrupole const & ellipse,
    afw::geom::Point2D const & center
) const {
    PTR(afw::detection::Footprint) region = afw::detection::growFootprint(
        footprint,
        getControl().region.nGrowFootprint,
        true
    );
    afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
    fullEllipse.getCore().scale(getControl().region.nInitialRadii);
    region = mergeFootprints(*region, afw::detection::Footprint(fullEllipse));
    double originalArea = region->getArea();
    region->clipTo(mask.getBBox(afw::image::PARENT));
    region->intersectMask(mask, _impl->badPixelMask);
    if ((1.0 - region->getArea() / originalArea) > getControl().region.maxBadPixelFraction) {
        region.reset();
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

    // Initialize the result object with NaNs and general failure flags set.
    Result result = _impl->makeResult();

    // Grow the footprint, clip bad pixels and the exposure bbox
    PTR(afw::detection::Footprint) initialFitRegion = determineInitialFitRegion(
        *exposure.getMaskedImage().getMask(),
        footprint
    );
    if (!initialFitRegion) {
        result.setFlag(CModelResult::MAX_BAD_PIXEL_FRACTION, true);
        return result;
    }
    if (initialFitRegion->getArea() > getControl().region.maxArea) {
        result.setFlag(CModelResult::MAX_AREA, true);
        return result;
    }

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
    _impl->initial.fit(getControl().initial, result.initial, initialData, exposure, *initialFitRegion);

    // Include a multiple of the initial-fit ellipse in the footprint, re-do clipping
    PTR(afw::detection::Footprint) finalFitRegion = determineFinalFitRegion(
        *exposure.getMaskedImage().getMask(),
        footprint,
        result.initial.ellipse,
        center
    );
    if (!finalFitRegion) {
        result.setFlag(CModelResult::MAX_BAD_PIXEL_FRACTION, true);
        return result;
    }
    if (finalFitRegion->getArea() > getControl().region.maxArea) {
        result.setFlag(CModelResult::MAX_AREA, true);
        return result;
    }

    // Do the exponential fit
    CModelStageData expData = initialData.changeModel(*_impl->exp.model);
    _impl->exp.fit(getControl().exp, result.exp, expData, exposure, *finalFitRegion);

    // Do the de Vaucouleur fit
    CModelStageData devData = initialData.changeModel(*_impl->dev.model);
    _impl->dev.fit(getControl().dev, result.dev, devData, exposure, *finalFitRegion);

    // Do the linear combination fit
    _impl->fitLinear(getControl(), result, expData, devData, exposure, *finalFitRegion);
    return result;
}

CModelAlgorithm::Result CModelAlgorithm::applyForced(
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    afw::geom::Point2D const & center,
    CModelResult const & reference,
    Scalar approxFlux
) const {

    // Negative approxFlux means we should come up with an estimate ourselves.
    // This is only used to avoid scaling problems in the optimizer, so it doesn't have to be very good.
    if (approxFlux < 0.0) {
        approxFlux = computeFluxInFootprint(*exposure.getMaskedImage().getImage(), footprint);
    }

    // Initialize the result object with NaNs and general failure flags set.
    Result result = _impl->makeResult();

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

    // Grow the footprint and include the initial ellipse, clip bad pixels and the exposure bbox;
    // in forced mode we can just use the final fit region immediately since we won't be changing
    // the initial fit ellipse.
    PTR(afw::detection::Footprint) finalFitRegion = determineFinalFitRegion(
        *exposure.getMaskedImage().getMask(),
        footprint,
        _impl->initial.ellipses.front().getCore(),
        center
    );
    if (!finalFitRegion) {
        result.setFlag(CModelResult::MAX_BAD_PIXEL_FRACTION, true);
        return result;
    }
    if (finalFitRegion->getArea() > getControl().region.maxArea) {
        result.setFlag(CModelResult::MAX_AREA, true);
        return result;
    }

    // Do the initial fit (amplitudes only)
    _impl->initial.fitLinear(getControl().initial, result.initial, initialData, exposure, *finalFitRegion);
    result.initial.setFlag(CModelStageResult::FAILED, false);

    // Do the exponential fit (amplitudes only)
    CModelStageData expData = initialData.changeModel(*_impl->exp.model);
    expData.nonlinear.deep() = reference.exp.nonlinear;
    expData.fixed.deep() = reference.exp.fixed;
    _impl->exp.fitLinear(getControl().exp, result.exp, expData, exposure, *finalFitRegion);

    // Do the de Vaucouleur fit (amplitudes only)
    CModelStageData devData = initialData.changeModel(*_impl->dev.model);
    devData.nonlinear.deep() = reference.dev.nonlinear;
    devData.fixed.deep() = reference.dev.fixed;
    _impl->dev.fitLinear(getControl().dev, result.dev, devData, exposure, *finalFitRegion);

    // Do the linear combination fit
    _impl->fitLinear(getControl(), result, expData, devData, exposure, *finalFitRegion);
    return result;
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
            pex::exceptions::LogicErrorException,
            "Algorithm was not initialized with a schema; cannot run in plugin mode"
        );
    }
    if (!exposure.getWcs()) {
        source.set(_impl->keys->flags[Result::NO_WCS], true);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Exposure has no Wcs"
        );
    }
    if (!exposure.getCalib()) {
        source.set(_impl->keys->flags[Result::NO_CALIB], true);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Exposure has no Calib"
        );
    }
    extensions::multiShapelet::FitPsfModel psfModel(*_impl->fitPsfCtrl, source);
    if (psfModel.hasFailed() || !(psfModel.ellipse.getArea() > 0.0)) {
        source.set(_impl->keys->flags[Result::NO_PSF], true);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Multishapelet PSF approximation failed or was not run"
        );
    }
    return psfModel.asMultiShapelet();
}

template <typename PixelT>
void CModelAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction psf = _processInputs(source, exposure);
    if (!source.getTable()->getShapeKey().isValid() ||
        (source.getTable()->getShapeFlagKey().isValid() && source.getShapeFlag())) {
        source.set(_impl->keys->flags[Result::NO_SHAPE], true);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Shape slot algorithm failed or was not run"
        );
    }
    // If PsfFlux has been run, use that for approx flux; otherwise we'll compute it ourselves.
    Scalar approxFlux = -1.0;
    if (source.getTable()->getPsfFluxKey().isValid() && !source.getPsfFluxFlag()) {
        approxFlux = source.getPsfFlux();
    }
    Result result = apply(exposure, *source.getFootprint(), psf, center, source.getShape(), approxFlux);
    _impl->keys->copyResultToRecord(result, source);
}

template <typename PixelT>
void CModelAlgorithm::_applyForced(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    afw::table::SourceRecord const & reference,
    afw::geom::AffineTransform const & refToMeas
) const {
    shapelet::MultiShapeletFunction psf = _processInputs(source, exposure);
    if (!_impl->refKeys) { // ideally we'd do this in the ctor, but we can't so we do it on first use
        _impl->refKeys.reset(
            new CModelKeys(
                *_impl->initial.model, *_impl->exp.model, *_impl->dev.model,
                reference.getSchema(), getControl().name
            )
        );
    }
    // If PsfFlux has been run, use that for approx flux; otherwise we'll compute it ourselves.
    Scalar approxFlux = -1.0;
    if (source.getTable()->getPsfFluxKey().isValid() && !source.getPsfFluxFlag()) {
        approxFlux = source.getPsfFlux();
    }
    Result refResult = _impl->refKeys->copyRecordToResult(source);
    Result result = applyForced(exposure, *source.getFootprint(), psf, center, refResult, approxFlux);
    _impl->keys->copyResultToRecord(result, source);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(CModelAlgorithm);

}}} // namespace lsst::meas::multifit