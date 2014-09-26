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

#ifndef LSST_MEAS_MULTIFIT_CModelFit_h_INCLUDED
#define LSST_MEAS_MULTIFIT_CModelFit_h_INCLUDED

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/algorithms/Algorithm.h"
#include "lsst/shapelet/RadialProfile.h"
#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/ProjectedLikelihood.h"
#include "lsst/meas/multifit/optimizer.h"

namespace lsst { namespace meas { namespace multifit {

class CModelAlgorithm;

struct CModelStageControl {

    CModelStageControl() :
        profileName("lux"),
        priorSource("CONFIG"),
        priorName(),
        nComponents(8),
        maxRadius(0),
        doRecordHistory(true),
        doRecordTime(true)
    {}

    shapelet::RadialProfile const & getProfile() const {
        return shapelet::RadialProfile::get(profileName);
    }

    PTR(Model) getModel() const;

    PTR(Prior) getPrior() const;

    LSST_CONTROL_FIELD(
        profileName, std::string,
        "Name of the shapelet.RadialProfile that defines the model to fit"
    );

    LSST_CONTROL_FIELD(
        priorSource, std::string,
        "One of 'FILE', 'CONFIG', or 'NONE', indicating whether the prior should be loaded from disk "
        "created from the nested prior config/control object, or None"
    );

    LSST_CONTROL_FIELD(
        priorName, std::string,
        "Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, "
        "with no extension), if priorSource='FILE'.  Ignored for forced fitting."
    );

    LSST_NESTED_CONTROL_FIELD(
        priorConfig, lsst.meas.multifit.multifitLib, SoftenedLinearPriorControl,
        "Configuration for the prior, used if priorSource='CONFIG'."
    );

    LSST_CONTROL_FIELD(nComponents, int, "Number of Gaussian used to approximate the profile");

    LSST_CONTROL_FIELD(
        maxRadius,
        int,
        "Maximum radius used in approximating profile with Gaussians (0=default for this profile)"
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.multifit.multifitLib, OptimizerControl,
        "Configuration for how the objective surface is explored.  Ignored for forced fitting"
    );

    LSST_NESTED_CONTROL_FIELD(
        likelihood, lsst.meas.multifit.multifitLib, ProjectedLikelihoodControl,
        "Configuration for how the compound model is evaluated and residuals are weighted in this "
        "stage of the fit"
    );

    LSST_CONTROL_FIELD(
        doRecordHistory, bool,
        "Whether to record the steps the optimizer takes (or just the number, if running as a plugin)"
    );

    LSST_CONTROL_FIELD(
        doRecordTime, bool,
        "Whether to record the time spent in this stage stage"
    );

};

struct CModelRegionControl {

    CModelRegionControl() :
        includePsfBBox(false),
        nGrowFootprint(5),
        nInitialRadii(3),
        maxArea(10000),
        maxBadPixelFraction(0.1)
    {
        badMaskPlanes.push_back("EDGE");
        badMaskPlanes.push_back("SAT");
    }

    LSST_CONTROL_FIELD(
        includePsfBBox, bool,
        "If True, always make the fit region at least the size of the PSF model realization's bounding box"
    );

    LSST_CONTROL_FIELD(
        nGrowFootprint, int,
        "Number of pixels to grow the original footprint by before the initial fit."
    );

    LSST_CONTROL_FIELD(
        nInitialRadii, double,
        "After the initial fit, extend the fit region to include all the pixels within "
        "this many initial-fit radii."
    );

    LSST_CONTROL_FIELD(
        maxArea, int,
        "Abort if the fit region grows beyond this many pixels."
    );

    LSST_CONTROL_FIELD(
        badMaskPlanes, std::vector<std::string>,
        "Mask planes that indicate pixels that should be ignored in the fit."
    );

    LSST_CONTROL_FIELD(
        maxBadPixelFraction, double,
        "Maximum fraction of pixels that may be ignored due to masks; "
        "more than this and we don't even try."
    );

};

struct CModelDiagnosticsControl {

    CModelDiagnosticsControl() : enabled(false), root(""), ids() {}

    LSST_CONTROL_FIELD(
        enabled,  bool,
        "Whether to write diagnostic outputs for post-run debugging"
    );

    LSST_CONTROL_FIELD(
        root, std::string,
        "Root output path for diagnostic outputs"
    );

    LSST_CONTROL_FIELD(
        ids, std::vector<boost::int64_t>,
        "Source IDs for which diagnostic outpust should be produced"
    );

};

struct CModelControl : public algorithms::AlgorithmControl {

    CModelControl() :
        algorithms::AlgorithmControl("cmodel", 2.5),
        psfName("multishapelet.psf"), minInitialRadius(0.1)
    {
        initial.nComponents = 3; // use very rough model in initial fit
        initial.optimizer.gradientThreshold = 1E-2; // with coarse convergence criteria
        initial.optimizer.minTrustRadiusThreshold = 1E-2;
        dev.profileName = "luv";
    }

    PTR(CModelControl) clone() const {
        return boost::static_pointer_cast<CModelControl>(_clone());
    }

    PTR(CModelAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        algorithms::AlgorithmMap const & others = algorithms::AlgorithmMap(),
        bool isForced = false
    ) const;

    LSST_CONTROL_FIELD(psfName, std::string, "Root name of the FitPsfAlgorithm fields.");

    LSST_NESTED_CONTROL_FIELD(
        region, lsst.meas.multifit.multifitLib, CModelRegionControl,
        "Configuration parameters related to the determination of the pixels to include in the fit."
    );

    LSST_NESTED_CONTROL_FIELD(
        diagnostics, lsst.meas.multifit.multifitLib, CModelDiagnosticsControl,
        "Configuration parameters related to diagnostic outputs for post-run debugging."
    );

    LSST_NESTED_CONTROL_FIELD(
        initial, lsst.meas.multifit.multifitLib, CModelStageControl,
        "An initial fit (usually with a fast, approximate model) used to warm-start the exp and dev fits, "
        "convolved with only the zeroth-order terms in the multi-shapelet PSF approximation."
    );

    LSST_NESTED_CONTROL_FIELD(
        exp, lsst.meas.multifit.multifitLib, CModelStageControl,
        "Independent fit of the exponential component"
    );

    LSST_NESTED_CONTROL_FIELD(
        dev, lsst.meas.multifit.multifitLib, CModelStageControl,
        "Independent fit of the de Vaucouleur component"
    );

    LSST_NESTED_CONTROL_FIELD(
        likelihood, lsst.meas.multifit.multifitLib, ProjectedLikelihoodControl,
        "configuration for how the compound model is evaluated and residuals are weighted in the exp+dev "
        "linear combination fit"
    );

    LSST_CONTROL_FIELD(
        minInitialRadius, double,
        "Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)"
    );

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata,
        algorithms::AlgorithmMap const & others,
        bool isForced
    ) const;
};

struct CModelStageResult {

    enum FlagBit {
        FAILED=0,
        TR_SMALL,
        MAX_ITERATIONS,
        NUMERIC_ERROR,
        N_FLAGS
    };

    CModelStageResult();

    PTR(Model) model;
    PTR(Prior) prior;
    PTR(OptimizerObjective) objfunc;
    Scalar flux;
    Scalar fluxSigma;
    Scalar objective;
    Scalar time;
    afw::geom::ellipses::Quadrupole ellipse;

    bool getFlag(FlagBit b) const { return flags[b]; }
    void setFlag(FlagBit b, bool value) { flags[b] = value; }

    ndarray::Array<Scalar const,1,1> nonlinear;
    ndarray::Array<Scalar const,1,1> amplitudes;
    ndarray::Array<Scalar const,1,1> fixed;

    afw::table::BaseCatalog history;
#ifndef SWIG
    std::bitset<N_FLAGS> flags;
#endif
};

struct CModelResult {

    enum FlagBit {
        FAILED=0,
        MAX_AREA,
        MAX_BAD_PIXEL_FRACTION,
        NO_SHAPE,
        NO_PSF,
        NO_WCS,
        NO_CALIB,
        N_FLAGS
    };

    CModelResult();

    Scalar flux;
    Scalar fluxSigma;
    Scalar fracDev;
    Scalar objective;

    bool getFlag(FlagBit b) const { return flags[b]; }
    void setFlag(FlagBit b, bool value) { flags[b] = value; }

    CModelStageResult initial;
    CModelStageResult exp;
    CModelStageResult dev;

    PTR(afw::detection::Footprint) initialFitRegion;
    PTR(afw::detection::Footprint) finalFitRegion;

#ifndef SWIG
    std::bitset<N_FLAGS> flags;
#endif
};

class CModelAlgorithm : public algorithms::Algorithm {
public:

    typedef CModelControl Control;
    typedef CModelResult Result;

    CModelAlgorithm(
        Control const & ctrl,
        afw::table::Schema & schema,
        algorithms::AlgorithmMap const & others,
        bool isForced
    );

    explicit CModelAlgorithm(Control const & ctrl);

    Control const & getControl() const {
        return static_cast<Control const &>(algorithms::Algorithm::getControl());
    }

    /**
     *  @brief Determine the initial fit region for a CModelAlgorithm fit
     *
     *  This routine grows the given footprint by nGrowFootprint, then clips on the bounding box
     *  of the given mask and removes pixels indicated as bad by badMaskPlanes.
     *
     *  If more than maxBadPixelFraction pixels are clipped, the returned pointer is null.
     */
    PTR(afw::detection::Footprint) determineInitialFitRegion(
        afw::image::Mask<> const & mask,
        afw::detection::Footprint const & footprint,
        afw::geom::Box2I const & psfBBox
    ) const;

    /**
     *  @brief Determine the final fit region for a CModelAlgorithm fit.
     *
     *  This routine grows the given footprint by nGrowFootprint, then extends it to include
     *  the given ellipse scaled by nInitialRadii.  It then clips on the bounding box of the
     *  given mask and removes pixels indicated as bad by badMaskPlanes.
     *
     *  If more than maxBadPixelFraction pixels are clipped, the returned pointer is null.
     */
    PTR(afw::detection::Footprint) determineFinalFitRegion(
        afw::image::Mask<> const & mask,
        afw::detection::Footprint const & footprint,
        afw::geom::Box2I const & psfBBox,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const;

    Result apply(
        afw::image::Exposure<Pixel> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        afw::geom::Point2D const & center,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar approxFlux=-1
    ) const;

    Result applyForced(
        afw::image::Exposure<Pixel> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        afw::geom::Point2D const & center,
        Result const & reference,
        Scalar approxFlux=-1
    ) const;

    void writeResultToRecord(Result const & result, afw::table::BaseRecord & record) const;

private:

    friend class CModelAlgorithmControl;

    // Actual implementations go here, so we can get partial results to the plugin
    // version when we throw.
    void _applyImpl(
        Result & result,
        afw::image::Exposure<Pixel> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        afw::geom::Point2D const & center,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar approxFlux
    ) const;

    // this method just throws an exception; it's present to resolve dispatches from _apply()
    void _applyImpl(
        Result & result,
        afw::image::Exposure<double> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        afw::geom::Point2D const & center,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar approxFlux=-1
    ) const {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "double-precision image measurement not implemented for CModelAlgorithm"
        );
    }

    // Actual implementations go here, so we can get partial results to the plugin
    // version when we throw.
    void _applyForcedImpl(
        Result & result,
        afw::image::Exposure<Pixel> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        afw::geom::Point2D const & center,
        Result const & reference,
        Scalar approxFlux
    ) const;

    // this method just throws an exception; it's present to resolve dispatches from _applyForced()
    void _applyForcedImpl(
        Result & result,
        afw::image::Exposure<double> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        afw::geom::Point2D const & center,
        Result const & reference,
        Scalar approxFlux=-1
    ) const {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "double-precision image measurement not implemented for CModelAlgorithm"
        );
    }

    template <typename PixelT>
    shapelet::MultiShapeletFunction _processInputs(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const;

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    template <typename PixelT>
    void _applyForced(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(CModelAlgorithm);

    class Impl;

    PTR(Impl) _impl;
};

inline PTR(CModelAlgorithm) CModelControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    return boost::static_pointer_cast<CModelAlgorithm>(
        _makeAlgorithm(schema, metadata, others, isForced)
    );
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_CModelFit_h_INCLUDED
