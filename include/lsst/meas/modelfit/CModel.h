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

#ifndef LSST_MEAS_MODELFIT_CModelFit_h_INCLUDED
#define LSST_MEAS_MODELFIT_CModelFit_h_INCLUDED

#include <bitset>
#include <string>

#include "ndarray.h"

#include "lsst/geom.h"
#include "lsst/pex/config.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/afw/table/Source.h"
#include "lsst/shapelet/RadialProfile.h"
#include "lsst/meas/modelfit/Model.h"
#include "lsst/meas/modelfit/Prior.h"
#include "lsst/meas/modelfit/MixturePrior.h"
#include "lsst/meas/modelfit/SoftenedLinearPrior.h"
#include "lsst/meas/modelfit/SemiEmpiricalPrior.h"
#include "lsst/meas/modelfit/UnitTransformedLikelihood.h"
#include "lsst/meas/modelfit/optimizer.h"
#include "lsst/meas/modelfit/PixelFitRegion.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  @page modelfitCModel CModel Magnitudes
 *
 *  The CModel approach to model-fit galaxy photometry - also known as the "Sloan Swindle" - is an
 *  approximation to bulge+disk or Sersic model fitting that follows the following sequence:
 *   - Fit a PSF-convolved elliptical exponential (Sersic n=1) model to the data.
 *   - Fit a PSF-convolved elliptical de Vaucouleurs (Sersic n=4) model to the data.
 *   - Holding the positions and ellipses of both models fixed (only allowing the amplitudes to vary),
 *     fit a linear combination of the two models.
 *  In the limit of pure bulge or pure disk galaxies, this approach yields the same results as a more
 *  principled bugle+disk or Sersic fit.  For galaxies that are a combination of the two components (or
 *  have more complicated morphologies, as of course all real galaxies do), it provides a smooth transition
 *  between the two models, and the fraction of flux in each of the two parameters is correlated with
 *  Sersic index and the true bulge-disk ratio.  Most importantly, this approach yieled good galaxy colors
 *  in the SDSS data processing.
 *
 *  In this implementation of the CModel algorithm, we actually have 4 stages:
 *   - In the "initial" stage, we fit a very approximate PSF-convolved elliptical model, just to provide
 *     a good starting point for the subsequence exponential and de Vaucouleur fits.  Because we use
 *     shapelet/Gaussian approximations to convolved models with the PSF, model evaluation is much faster
 *     when only a few Gaussians are used in the approximation, as is done here.  In the future, we may
 *     also use a simpler PSF approximation in the initial fit, but this is not yet implemented.  We also
 *     have not yet researched how best to make use of the initial fit (i.e. how does the initial best-fit
 *     radius typically relate to the best-fit exponential radius?), or what convergence criteria should
 *     be used in the initial fit.  Following the initial fit, we also revisit the question of which pixels
 *     should be included in the fit (see CModelRegionControl).
 *   - In the "exp" stage, we start with the "initial" fit results, and fit an elliptical exponential
 *     profile.
 *   - In the "dev" stage, we start with the "initial" fit results, and fit an elliptical de Vaucouleur
 *     profile.
 *   - Holding the "exp" and "dev" ellipses fixed, we fit a linear combination of those two profiles.
 *  In all of these steps, the centroid is held fixed at a given input value (take from the slot centroid
 *  when run by the measurement framework).
 *
 *  @section cmodelUnits Units
 *
 *  Unlike most measurement algorithms, CModel requires the Exposure it is given to have both a Wcs and
 *  a PhotoCalib.  This is because it makes use of Bayesian priors, and hence it has to know the relationship
 *  between the raw units of the image (pixel and count) and the global units in which the priors are defined.
 *
 *  In fact, all of the nonlinear fits in CModel are done in a special, local coordinate system, defined
 *  by a Wcs in which the "pixels" have units of arcseconds (because we never create an image in this system,
 *  we don't have to worry about the size of the pixels) and the fluxes should be of order unity.  In
 *  addition to allowing us to use priors, it also ensures that the parameters all have the same order
 *  of magnitude, which improves the behavior of the optimizer.
 *
 *  See @ref modelfitUnits for more information.
 *
 *  @section cmodelForced Forced Photometry
 *
 *  In forced photometry, we replace the three nonlinear fits with amplitude-only fits, and then repeat the
 *  final linear fit, using the ellipses from the reference catalog in all casees.  We do allow the relative
 *  amplitudes of the two components to vary in forced mode, though in the future we will add an option to
 *  hold this fixed as well as the ellipses.
 *
 *  @section cmodelPsf Shapelet Approximations to the PSF
 *
 *  The CModel algorithm relies on a multi-shapelet approximation to the PSF to convolve galaxy models.  It
 *  does not compute this approximation directly; for CModelAlgorithm methods that take inputs directly
 *  as arguments, the PSF must be supplied as a shapelet::MultiShapeletFunction instance.  When using
 *  SourceRecords for input/output, CModel assumes that the ShapeletPsfApprox plugin has already been
 *  run (see psf.py), and uses the fields created by that plugin to retrieve the PSF approximation.
 *
 *  @section cmodelOrg Code Organization
 *
 *  The CModel implementation consists of many classes, defined in this file and CModel.cc.  These mostly
 *  fall into four categories:
 *    - Control structs: C++ analogs of Python Config classes, these structs contain the
 *      configuration parameters that control the behavior of the algorithm.  These are nested; the
 *      @ref CModelControl struct contains a @ref PixelFitRegionControl instance and three
 *      @ref CModelStageControl (one for each of "initial", "exp", and "dev").  The configuration
 *      for the final amplitude-only fit goes in @ref CModelControl itself; because it is a simpler
 *      linear fit, it doesn't have much in common with the first three stages.
 *    - Result structs: while the algorithm has methods to use SourceRecord objects for input/output,
 *      it can also take inputs directly as arguments and return the outputs using these structs.  Like
 *      the Control structs, the master @ref CModelResult struct holds three @ref CModelStageResult classes,
 *      for each of the three nonlinear fits.
 *    - Keys structs: these private classes (defined in an anonymous namespace in CModel.cc) hold the
 *      afw::table::Key and FunctorKey objects that provide a mapping from the Result structs to
 *      Schema fields.  They also provide methods to transfer values from Results to Records, or the
 *      reverse.  These are also split into a master CModelKeys struct, which holds three
 *      CModelStageKeys structs.
 *    - Impl classes: these private classes contain the actual algorithmic code.  Once again, we
 *      have the master implementation class (CModelAlgorithm::Impl) and a class for the nonlinear
 *      fitting stages (CModelStageImpl).
 *  In addition to these categories, we also have the @ref CModelAlgorithm class, which is the C++ public
 *  interface to all of this, and the CModelStageData class, a private class that aggregates per-source
 *  state and makes it easier to pass it around.
 */

/**
 *  Nested control object for CModel that configures one of the three ("initial", "exp", "dev") nonlinear
 *  fitting stages.
 */
struct CModelStageControl {

    CModelStageControl() :
        profileName("lux"),
        priorSource("EMPIRICAL"),
        priorName(),
        nComponents(8),
        maxRadius(0),
        usePixelWeights(false),
        weightsMultiplier(1.0),
        doRecordHistory(true),
        doRecordTime(true)
    {}

    shapelet::RadialProfile const & getProfile() const {
        return shapelet::RadialProfile::get(profileName);
    }

    std::shared_ptr<Model> getModel() const;

    std::shared_ptr<Prior> getPrior() const;

    LSST_CONTROL_FIELD(
        profileName, std::string,
        "Name of the shapelet.RadialProfile that defines the model to fit"
    );

    LSST_CONTROL_FIELD(
        priorSource, std::string,
        "One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded "
        "from disk, created from one of the nested prior config/control objects, or None"
    );

    LSST_CONTROL_FIELD(
        priorName, std::string,
        "Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, "
        "with no extension), if priorSource='FILE'.  Ignored for forced fitting."
    );

    LSST_NESTED_CONTROL_FIELD(
        linearPriorConfig, lsst.meas.modelfit.priors, SoftenedLinearPriorControl,
        "Configuration for a linear prior, used if priorSource='LINEAR'."
    );

    LSST_NESTED_CONTROL_FIELD(
        empiricalPriorConfig, lsst.meas.modelfit.priors, SemiEmpiricalPriorControl,
        "Configuration for an empirical prior, used if priorSource='EMPIRICAL'."
    );

    LSST_CONTROL_FIELD(nComponents, int, "Number of Gaussian used to approximate the profile");

    LSST_CONTROL_FIELD(
        maxRadius,
        int,
        "Maximum radius used in approximating profile with Gaussians (0=default for this profile)"
    );

    LSST_CONTROL_FIELD(
        usePixelWeights,
        bool,
        "Use per-pixel variances as weights in the nonlinear fit (the final linear fit for"
        " flux never uses per-pixel variances)"
    );

    LSST_CONTROL_FIELD(
        weightsMultiplier,
        double,
        "Scale the likelihood by this factor to artificially reweight it w.r.t. the prior."
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.modelfit.optimizer, OptimizerControl,
        "Configuration for how the objective surface is explored.  Ignored for forced fitting"
    );

    LSST_CONTROL_FIELD(
        doRecordHistory, bool,
        "Whether to record the steps the optimizer takes (or just the number, if running as a plugin)"
    );

    LSST_CONTROL_FIELD(
        doRecordTime, bool,
        "Whether to record the time spent in this stage"
    );

};

/**
 *  The main control object for CModel, containing parameters for the final linear fit and aggregating
 *  the other control objects.
 */
struct CModelControl {

    CModelControl() :
        psfName("modelfit_DoubleShapeletPsfApprox"),
        minInitialRadius(0.1),
        fallbackInitialMomentsPsfFactor(1.5)
    {
        initial.nComponents = 3; // use very rough model in initial fit
        initial.optimizer.gradientThreshold = 1E-3; // with slightly coarser convergence criteria
        initial.optimizer.minTrustRadiusThreshold = 1E-2;
        initial.usePixelWeights = true;
        dev.profileName = "luv";
        exp.nComponents = 6;
        exp.optimizer.maxOuterIterations = 250;
    }

    LSST_CONTROL_FIELD(
        psfName,
        std::string,
        "Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; "
        "must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey."
    );

    LSST_NESTED_CONTROL_FIELD(
        region, lsst.meas.modelfit.pixelFitRegion, PixelFitRegionControl,
        "Configuration parameters related to the determination of the pixels to include in the fit."
    );

    LSST_NESTED_CONTROL_FIELD(
        initial, lsst.meas.modelfit.cmodel, CModelStageControl,
        "An initial fit (usually with a fast, approximate model) used to warm-start the exp and dev fits, "
        "convolved with only the zeroth-order terms in the multi-shapelet PSF approximation."
    );

    LSST_NESTED_CONTROL_FIELD(
        exp, lsst.meas.modelfit.cmodel, CModelStageControl,
        "Independent fit of the exponential component"
    );

    LSST_NESTED_CONTROL_FIELD(
        dev, lsst.meas.modelfit.cmodel, CModelStageControl,
        "Independent fit of the de Vaucouleur component"
    );

    LSST_CONTROL_FIELD(
        minInitialRadius, double,
        "Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)"
    );

    LSST_CONTROL_FIELD(
        fallbackInitialMomentsPsfFactor, double,
        "If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this."
        "  If <= 0.0, abort the fit early instead."
    );

};

/**
 *  Result object for a single nonlinear fitting stage of the CModel algorithm
 */
struct CModelStageResult {

    /// Flags for a single CModel stage (note that there are additional flags for the full multi-stage fit)
    enum FlagBit {
        FAILED=0,        ///< General flag, indicating whether the flux for this stage can be trusted.
        TR_SMALL,        ///< Whether convergence was due to the optimizer trust region getting too small
                         ///  (not a failure!)
        MAX_ITERATIONS,  ///< Whether the optimizer exceeded the maximum number of iterations.  Indicates
                         ///  a suspect fit, but not necessarily a bad one (implies FAILED).
        NUMERIC_ERROR,   ///< Optimizer encountered a numerical error (something likely went to infinity).
                         ///  Result will be unusable; implies FAILED.
        BAD_REFERENCE,   ///< Reference fit failed, so forced fit will fail as well.
        NO_FLUX,         ///< No flux was measured.
        N_FLAGS          ///< Non-flag counter to indicate the number of flags
    };

    CModelStageResult();

    std::shared_ptr<Model> model;    ///< Model object that defines the parametrization (defined fully by Control struct)
    std::shared_ptr<Prior> prior;    ///< Bayesian priors on the parameters (defined fully by Control struct)
    std::shared_ptr<OptimizerObjective> objfunc;  ///< Objective class used by the optimizer
    std::shared_ptr<UnitTransformedLikelihood> likelihood; ///< Object used to evaluate models and compare to data.
    Scalar instFlux;         ///< Flux measured from just this stage fit.
    Scalar instFluxErr;    ///< Flux uncertainty from just this stage fit.
    Scalar instFluxInner;    ///< Flux measured strictly within the fit region (no extrapolation).
    Scalar objective;    ///< Value of the objective function at the best fit point: chisq/2 - ln(prior)
    Scalar time;         ///< Time spent in this fit in seconds.
    afw::geom::ellipses::Quadrupole ellipse;  ///< Best fit half-light ellipse in pixel coordinates

    ndarray::Array<Scalar const,1,1> nonlinear;  ///< Opaque nonlinear parameters in specialized units
    ndarray::Array<Scalar const,1,1> amplitudes; ///< Opaque linear parameters in specialized units
    ndarray::Array<Scalar const,1,1> fixed;      ///< Opaque fixed parameters in specialized units

    afw::table::BaseCatalog history;  ///< Trace of the optimizer's path, if enabled by diagnostic options
    std::bitset<N_FLAGS> flags; ///< Array of flags.
};

/**
 *  Master result object for CModel, containing results for the final linear fit and three nested
 *  CModelStageResult objects for the results of the previous stages.
 */
struct CModelResult {

    /// Flags that apply to all four CModel fits or just the last one.
    enum FlagBit {
        FAILED=0,                ///< General failure flag for the linear fit flux; set if any other
                                 ///  CModel flag is set, or if any of the three previous stages failed.
        REGION_MAX_AREA,                  ///< Set if we aborted early because the fit region was too large.
        REGION_MAX_BAD_PIXEL_FRACTION,    ///< Set if we aborted early because the fit region had too many
                                          ///  bad pixels.
        REGION_USED_FOOTPRINT_AREA,       ///< Kron radius was unavailable or outside bounds, so the
                                          ///  second-moment ellipse scaled to the footprint area was used
                                          ///  instead.
        REGION_USED_PSF_AREA,             ///< Kron radius was unavailable or outside bounds, so the
                                          ///  second-moment ellipse scaled to the PSF area was used instead.
        REGION_USED_INITIAL_ELLIPSE_MIN,  ///< Fit region implied by the best-fit ellipse of the initial was
                                          ///  too small, so we used the configuration minimum instead.
        REGION_USED_INITIAL_ELLIPSE_MAX,  ///< Fit region implied by the best-fit ellipse of the initial was
                                          ///  too large, so we used the configuration maximum instead.
        NO_SHAPE,                ///< Set if the input SourceRecord had no valid shape slot with which to
                                 ///  start the fit.
        SMALL_SHAPE,             ///< Initial moments were sufficiently small that we used minInitialRadius
                                 ///  to set the initial parameters.
        NO_SHAPELET_PSF,         ///< Set if the Psf shapelet approximation failed.
        BAD_CENTROID,            ///< Input centroid did not land within the fit region.
        BAD_REFERENCE,           ///< Reference fit failed, so forced fit will fail as well.
        NO_FLUX,                 ///< No flux was measured.
        N_FLAGS                  ///< Non-flag counter to indicate the number of flags
    };

    CModelResult();

    Scalar instFlux;       ///< Flux from the final linear fit
    Scalar instFluxErr;  ///< Flux uncertainty from the final linear fit
    Scalar instFluxInner;  ///< Flux measured strictly within the fit region (no extrapolation).
    Scalar fracDev;    ///< Fraction of flux from the final linear fit in the de Vaucouleur component
                       ///  (always between 0 and 1).
    Scalar objective;  ///< Objective value at the best-fit point (chisq/2)

    CModelStageResult initial; ///< Results from the initial approximate nonlinear fit that feeds the others
    CModelStageResult exp;     ///< Results from the exponential (Sersic n=1) fit
    CModelStageResult dev;     ///< Results from the de Vaucouleur (Sersic n=4) fit

    afw::geom::ellipses::Quadrupole initialFitRegion;  ///< Pixels used in the initial fit.
    afw::geom::ellipses::Quadrupole finalFitRegion;    ///< Pixels used in the exp, dev, and linear fits.

    LocalUnitTransform fitSysToMeasSys; ///< Transforms to the coordinate system where parameters are defined
    std::bitset<N_FLAGS> flags; ///< Array of flags.
};

/**
 *  Main public interface class for CModel algorithm.
 *
 *  See @ref modelfitCModel for a full description of the algorithm.
 *
 *  This class provides the methods that actually execute the algorithm, and (depending on how it is
 *  constructed) holds the Key objects necessary to use SourceRecords for input and output.
 */
class CModelAlgorithm {
public:

    typedef CModelControl Control; ///< Typedef to the master Control struct
    typedef CModelResult Result;   ///< Typedef to the master Result struct

    /**
     *  Construct an algorithm instance and add its fields to the Schema.
     *
     *  All fields needed to write the outputs of a regular, non-forced fit will be added to the given
     *  Schema.  In addition, keys needed to retrieve the PSF shapelet approximation (assuming the
     *  ShapeletPsfApprox plugin has been run) will be extracted from the Schema.
     *
     *  @param[in]     name    Name of the algorithm used as a prefix for all fields added to the Schema.
     *  @param[in]     ctrl    Control object that configures the algorithm.
     *  @param[in,out] schema  Schema to which fields will be added, and from which keys for the PSF
     *                         shapelet approximation will be extacted.
     */
    CModelAlgorithm(
        std::string const & name,
        Control const & ctrl,
        afw::table::Schema & schema
    );

    /**
     *  Construct an algorithm instance suitable for forced photometry and add its fields to the Schema.
     *
     *  All fields needed to write the outputs of a forced fit will be added to the given SchemaMapper's
     *  output schema.  Keys needed to retrieve the reference ellipses for the exp and dev fits will be
     *  extracted from the SchemaMapper's input schema.  In addition, keys needed to retrieve the PSF
     *  shapelet approximation (assuming the ShapeletPsfApprox plugin has been run) will be extracted
     *  from the SchemaMapper's output schema (note that the ShapeletPsfApprox plugin must be run in
     *  forced mode as well, to approximate the measurement image's PSF rather than the reference image's
     *  PSF, so its outputs are found in the output schema, not the input schema).
     *
     *  @param[in]     name    Name of the algorithm used as a prefix for all fields added to the Schema.
     *  @param[in]     ctrl    Control object that configures the algorithm.
     *  @param[in,out] schemaMapper  SchemaMapper containing input (reference) and output schemas.
     */
    CModelAlgorithm(
        std::string const & name,
        Control const & ctrl,
        afw::table::SchemaMapper & schemaMapper
    );

    /**
     *  Construct an algorithm instance that cannot use SourceRecords for input/output.
     *
     *  This constructor initializes the algorithm without initializing any of the keys necessary to
     *  operate on SourceRecords.  As a result, only methods that take inputs directly and return Result
     *  objects may be called.
     */
    explicit CModelAlgorithm(Control const & ctrl);

    /// Return the control object the algorithm was constructed with.
    Control const & getControl() const { return _ctrl; }

    /**
     *  Run the CModel algorithm on an image, supplying inputs directly and returning outputs in a Result.
     *
     *  @param[in]   exposure     Image to measure.  Must have a valid Psf, Wcs and PhotoCalib.
     *  @param[in]   psf          multi-shapelet approximation to the PSF at the position of the source
     *  @param[in]   center       Centroid of the source to be fit.
     *  @param[in]   moments      Non-PSF-corrected moments of the source, used to initialize the model
     *                            parameters
     *  @param[in]   approxFlux   Rough estimate of the flux of the source, used to set the fit coordinate
     *                            system and ensure internal parameters are of order unity.  If less than
     *                            or equal to zero, the sum of the flux within the footprint will be used.
     *  @param[in]   kronRadius   Estimate of the Kron radius (optional); used as the first choice when
     *                            estimating the region of pixel to include in the fit.
     *  @param[in]   footprintArea  Area of the detection Fooptrint; used as the fallback when
     *                              estimating the region of pixel to include in the fit.
     */
    Result apply(
        afw::image::Exposure<Pixel> const & exposure,
        shapelet::MultiShapeletFunction const & psf,
        geom::Point2D const & center,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar approxFlux=-1,
        Scalar kronRadius=-1,
        int footprintArea=-1
    ) const;

    /**
     *  Run the CModel algorithm in forced mode on an image, supplying inputs directly and returning
     *  outputs in a Result.
     *
     *  @param[in]   exposure     Image to measure.  Must have a valid Psf, Wcs and PhotoCalib.
     *  @param[in]   psf          multi-shapelet approximation to the PSF at the position of the source
     *  @param[in]   center       Centroid of the source to be fit.
     *  @param[in]   reference    Result object from a previous, non-forced run of CModelAlgorithm.
     *  @param[in]   approxFlux   Rough estimate of the flux of the source, used to set the fit coordinate
     *                            system and ensure internal parameters are of order unity.  If less than
     *                            or equal to zero, the sum of the flux within the footprint will be used.
     */
    Result applyForced(
        afw::image::Exposure<Pixel> const & exposure,
        shapelet::MultiShapeletFunction const & psf,
        geom::Point2D const & center,
        Result const & reference,
        Scalar approxFlux=-1
    ) const;

    /**
     *  Run the CModel algorithm on an image, using a SourceRecord for inputs and outputs.
     *
     *  @param[in,out] measRecord  A SourceRecord instance used to provide a Footprint, the centroid and
     *                             shape of the source, a MultiShapeletFunction PSF, and an approximate
     *                             estimate of the (via the PsfFlux slot), and to which all outputs will
     *                             be written.
     *  @param[in]     exposure    Image to be measured.  Must have a valid Psf, Wcs, and PhotoCalib.
     *
     *  To run this method, the CModelAlgorithm instance must have been created using the constructor
     *  that takes a Schema argument, and that Schema must match the Schema of the SourceRecord passed here.
     */
     void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<Pixel> const & exposure
    ) const;

    /**
     *  Run the CModel algorithm in forced mode on an image, using a SourceRecord for inputs and outputs.
     *
     *  @param[in,out] measRecord  A SourceRecord instance used to provide a Footprint, the centroid of
     *                             the source, a MultiShapeletFunction PSF, and an approximate
     *                             estimate of the (via the PsfFlux slot), and to which all outputs will
     *                             be written.
     *  @param[in]     exposure    Image to be measured.  Must have a valid Psf, Wcs, and PhotoCalib.
     *  @param[in]     refRecord   A SourceRecord that contains the outputs of a previous non-forced run
     *                             of CModelAlgorithm (which may have taken place on an image with a
     *                             different Wcs).
     *
     *  To run this method, the CModelAlgorithm instance must have been created using the constructor
     *  that takes a Schema argument, and that Schema must match the Schema of the SourceRecord passed here.
     */
    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<Pixel> const & exposure,
        afw::table::SourceRecord const & refRecord
    ) const;

    /**
     *  Handle an exception thrown by one of the measure() methods, setting the appropriate flag in
     *  the given record.
     *
     *  @param[out]   measRecord   Record on which the flag should be set.
     *  @param[in]    error        Error containing the bit to be set.  If null, only the general
     *                             failure bit will be set.
     */
    void fail(
        afw::table::SourceRecord & measRecord,
        meas::base::MeasurementError * error
    ) const;

    /// Copy values from a Result struct to a BaseRecord object.
    void writeResultToRecord(Result const & result, afw::table::BaseRecord & record) const;

private:

    friend class CModelAlgorithmControl;

    // Actual implementations go here; we use an output argument for the result so we can get partial
    // results to the plugin version when we throw.
    void _applyImpl(
        Result & result,
        afw::image::Exposure<Pixel> const & exposure,
        shapelet::MultiShapeletFunction const & psf,
        geom::Point2D const & center,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar approxFlux,
        Scalar kronRadius=-1,
        int footprintArea=-1
    ) const;

    // Actual implementations go here; we use an output argument for the result so we can get partial
    // results to the SourceRecord version when we throw.
    void _applyForcedImpl(
        Result & result,
        afw::image::Exposure<Pixel> const & exposure,
        shapelet::MultiShapeletFunction const & psf,
        geom::Point2D const & center,
        Result const & reference,
        Scalar approxFlux
    ) const;

    // gets/checks inputs from SourceRecord that are needed by both apply and applyForced
    template <typename PixelT>
    shapelet::MultiShapeletFunction _processInputs(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure
    ) const;

    class Impl;

    Control _ctrl;
    std::shared_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_CModelFit_h_INCLUDED
