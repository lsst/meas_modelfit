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

#ifndef LSST_MEAS_MODELFIT_GeneralPsfFitter_h_INCLUDED
#define LSST_MEAS_MODELFIT_GeneralPsfFitter_h_INCLUDED

#include <memory>

#include "lsst/pex/config.h"
#include "lsst/shapelet/FunctorKeys.h"
#include "lsst/meas/modelfit/Model.h"
#include "lsst/meas/modelfit/Prior.h"
#include "lsst/meas/modelfit/Likelihood.h"
#include "lsst/geom.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/modelfit/optimizer.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  Control object used to define one piece of multishapelet fit to a PSF model; see GeneralPsfFitterControl
 */
class GeneralPsfFitterComponentControl {
public:

    GeneralPsfFitterComponentControl(int order_=0, double radiusFactor_=1.0) :
        order(order_), positionPriorSigma(0.1), ellipticityPriorSigma(0.3),
        radiusFactor(radiusFactor_), radiusPriorSigma(0.5)
    {}

    LSST_CONTROL_FIELD(
        order, int,
        "shapelet order for this component; negative to disable this component completely"
    );
    LSST_CONTROL_FIELD(
        positionPriorSigma, double,
        "sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, "
        "relative to the center of the PSF image"
    );
    LSST_CONTROL_FIELD(
        ellipticityPriorSigma, double,
        "sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta"
    );
    LSST_CONTROL_FIELD(
        radiusFactor, double,
        "Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either "
        "the second-moments radius of the PSF image (in an initial fit), or the radius of the primary "
        "component in a previous fit.  Ignored if the previous fit included this component (as then we "
        "can just use that radius)."
    );
    LSST_CONTROL_FIELD(
        radiusPriorSigma, double,
        "sigma in a Gaussian prior on ln(radius/fiducialRadius)"
    );

};

/**
 * Control object used to configure a multishapelet fit to a PSF model; see GeneralPsfFitter.
 *
 *  The default configuration corresponds to fitting an elliptical double-Gaussian, in which each component
 *  can have different radii, positions, and ellipticities.  While the fitter can support much more complex
 *  models, at present, fitting these is prohibitively slow, and is not recommended in production
 *  environments (use DoubleShapeletPsfApprox instead).
 */
class GeneralPsfFitterControl {
public:

    GeneralPsfFitterControl() :
        inner(-1, 0.5), primary(0, 1.0), wings(0, 2.0), outer(-1, 4.0), defaultNoiseSigma(0.001)
    {}

    LSST_NESTED_CONTROL_FIELD(
        inner, lsst.meas.modelfit.psf.psfContinued, GeneralPsfFitterComponentControl,
        "Innermost shapelet expansion, used to fit PSFs with very sharp cores"
    );

    LSST_NESTED_CONTROL_FIELD(
        primary, lsst.meas.modelfit.psf.psfContinued, GeneralPsfFitterComponentControl,
        "Primary shapelet expansion, typically used to fit the bulk of the PSF "
    );

    LSST_NESTED_CONTROL_FIELD(
        wings, lsst.meas.modelfit.psf.psfContinued, GeneralPsfFitterComponentControl,
        "Wing shapelet expansion (between primary and outer), typically used to fit the wings of the PSF"
    );

    LSST_NESTED_CONTROL_FIELD(
        outer, lsst.meas.modelfit.psf.psfContinued, GeneralPsfFitterComponentControl,
        "Outermost shapelet expansion, used to fit PSFs with very broad wings"
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.modelfit.optimizer, OptimizerControl,
        "Configuration of the optimizer used to do the fitting"
    );

    LSST_CONTROL_FIELD(
        defaultNoiseSigma, double, "Default value for the noiseSigma parameter in GeneralPsfFitter.apply()"
    );

};

/**
 *  @brief Class for fitting multishapelet models to PSF images
 *
 *  This class fits up to four shapelet expansions simultaneously to a PSF image, with the relative radii
 *  and number of shapelet coefficients for each expansion separately configurable.  These expansions are
 *  also named; this allows us to map different fits with some expansions disabled to each other, in order
 *  to first fit an approximate model and follow this up with a more complete model, using the approximate
 *  model as a starting point.
 *
 *  The configuration also defines a simple Bayesian prior for the fit, defined using simple independent
 *  Gaussians for the ellipse parameters of each component.  The priors can be disabled by setting their
 *  width (xxPriorSigma in the control object) to infinity, and those parameters can be held fixed at
 *  their input values by setting the prior width to zero.  The priors are always centered at the input
 *  value, meaning that it may be more appropriate to think of the priors as a form of regularization,
 *  rather than a rigorous prior.  In fact, it's impossible to use a prior here rigorously without a
 *  noise model for the PSF image, which is something the LSST Psf class doesn't provide, and here is
 *  just provided as a constant noise sigma to be provided by the user (who generally just has to chose
 *  a small number arbitrarily).  Decreasing the noise sigma will of course decrease the effect of the
 *  priors (and vice versa).  In any case, having some sort of regularization is probably a good idea,
 *  as this is a very high-dimensional fit.
 */
class GeneralPsfFitter {
public:
    /// Initialize the fitter class with the given control object.
    explicit GeneralPsfFitter(GeneralPsfFitterControl const & ctrl);

    /**
     *  Add fields to a Schema that can be used to store the MultiShapeletFunction returned by apply().
     *
     *  @param[in,out]   schema    Schema to add fields to.
     *  @param[in]       prefix    Field name prefix for all fields.
     *  @return a FunctorKey that can get/set MultiShapeletFunctions that match the configuration of this
     *          fitter on a record.
     */
    shapelet::MultiShapeletFunctionKey addFields(
        afw::table::Schema & schema,
        std::string const & prefix
    ) const;

    /**
     *  Return the Model object that corresponds to the configuration.
     *
     *  In addition to the shapelet coefficients (stored in the "amplitudes" array), this Model
     *  stores all the initial ellipse parameters in the "fixed" array, as these are used to
     *  define the center of the prior; the "nonlinear" parameters are the free-to-vary ellipse
     *  parameters minus the corresponding initial values.
     */
    std::shared_ptr<Model> getModel() const { return _model; }

    /**
     *  Return the Prior object that corresponds to the configuration.
     *
     *  This Prior class only supports evaluate() and evaluateDerivatives(), reflecting the fact
     *  that we only intend to use it with a Optimizer, not a Sampler.
     */
    std::shared_ptr<Prior> getPrior() const { return _prior; }

    /**
     *  Adapt a differently-configured previous fit to be used as an starting point for this GeneralPsfFitter.
     *
     *  @param[in] previousFit     The return value of apply() from a differently-configured
     *                             instance of GeneralPsfFitter.
     *  @param[in] previousModel   The Model associated with the GeneralPsfFitter used to create previousFit.
     *
     *  @return a new MultiShapelet function that may be passed directly to apply().  When possible,
     *  the ellipse and shapelet coefficeints will be copied from previousFit; higher-order coefficients
     *  will be set to zero, and any components used in this but unused in the previous fit will have their
     *  ellipses set relative to the previous fit's "primary" component.
     */
    shapelet::MultiShapeletFunction adapt(
        shapelet::MultiShapeletFunction const & previousFit,
        std::shared_ptr<Model> previousModel
    ) const;

    //@{
    /**
     *  Perform an initial fit to a PSF image.
     *
     *  @param[in]  image       The image to fit, typically the result of Psf::computeKernelImage().  The
     *                          image's xy0 should be set such that the center of the PSF is at (0,0).
     *  @param[in]  moments     Second moments of the PSF, typically result of Psf::computeShape() or running
     *                          some other adaptive moments code on the PSF image.  This will be used to
     *                          set the initial ellipses of the multishapelet model.
     *  @param[in]  noiseSigma  An estimate of the noise in the image.  As LSST PSF images are generally
     *                          assumed to be noise-free, this is really just a fiddle-factor for the user.
     *                          A default value from the control object is used if this is negative.
     *  @param[in]  pState      Pointer to an integer which is used to return the optimizerState from apply.
     */
    shapelet::MultiShapeletFunction apply(
        afw::image::Image<Pixel> const & image,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar noiseSigma=-1,
        int * pState = nullptr
    ) const;
    shapelet::MultiShapeletFunction apply (
        afw::image::Image<double> const & image,
        afw::geom::ellipses::Quadrupole const & moments,
        Scalar noiseSigma=-1,
        int * pState = nullptr
    ) const {
        return apply(afw::image::Image<float>(image, true), moments, noiseSigma, pState);
    }
    //@}

    //@{
    /**
     *  Perform a fit to a PSF image, using a previous fit as a starting point
     *
     *  @param[in]  image       The image to fit, typically the result of Psf::computeKernelImage().  The
     *                          image's xy0 should be set such that the center of the PSF is at (0,0).
     *  @param[in]  initial     The result of a previous call to apply(), using an identically-configured
     *                          GeneralPsfFitter instance.  To use a result from a differently-configured GeneralPsfFitter,
     *                          use adapt().
     *  @param[in]  noiseSigma  An estimate of the noise in the image.  As LSST PSF images are generally
     *                          assumed to be noise-free, this is really just a fiddle-factor for the user.
     *                          A default value from the control object is used if this is negative.
     *  @param[in]  pState      Pointer to an integer which is used to return the optimizerState from apply.
     */
    shapelet::MultiShapeletFunction apply(
        afw::image::Image<Pixel> const & image,
        shapelet::MultiShapeletFunction const & initial,
        Scalar noiseSigma=-1,
        int * pState = nullptr
    ) const;
    shapelet::MultiShapeletFunction apply(
        afw::image::Image<double> const & image,
        shapelet::MultiShapeletFunction const & initial,
        Scalar noiseSigma=-1,
        int * pState = nullptr
    ) const {
        return apply(afw::image::Image<float>(image, true), initial, noiseSigma);
    }
    //@}

private:
    GeneralPsfFitterControl _ctrl;
    std::shared_ptr<Model> _model;
    std::shared_ptr<Prior> _prior;
};

class GeneralPsfFitterAlgorithm : public GeneralPsfFitter {
public:

    // Structures and routines to manage flaghandler
    static base::FlagDefinitionList const & getFlagDefinitions();
    static base::FlagDefinition const FAILURE;
    static base::FlagDefinition const MAX_INNER_ITERATIONS;
    static base::FlagDefinition const MAX_OUTER_ITERATIONS;
    static base::FlagDefinition const EXCEPTION;
    static base::FlagDefinition const CONTAINS_NAN;

    typedef GeneralPsfFitterControl Control;

    GeneralPsfFitterAlgorithm(GeneralPsfFitterControl const & ctrl,
        afw::table::Schema & schema,
        std::string const & prefix
    );

    shapelet::MultiShapeletFunctionKey getKey() {
        return _key;
    }

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Image<double> const & image,
        shapelet::MultiShapeletFunction const & initial
    ) const;

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Image<double> const & image,
        afw::geom::ellipses::Quadrupole const & moments
    ) const;

    void fail(
        afw::table::SourceRecord & measRecord,
        lsst::meas::base::MeasurementError * error=nullptr
    ) const;

private:
    shapelet::MultiShapeletFunctionKey _key;
    lsst::meas::base::FlagHandler _flagHandler;
};

/**
 *  Likelihood object used to fit multishapelet models to PSF model images; mostly for internal use
 *  by GeneralPsfFitter.
 */
class MultiShapeletPsfLikelihood : public Likelihood {
public:

    MultiShapeletPsfLikelihood(
        ndarray::Array<Pixel const,2,1> const & image,
        geom::Point2I const & xy0,
        std::shared_ptr<Model> model,
        Scalar sigma,
        ndarray::Array<Scalar const,1,1> const & fixed
    );

    void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        bool doApplyWeights=true
    ) const override;

    virtual ~MultiShapeletPsfLikelihood();

private:
    class Impl;
    std::unique_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_GeneralPsfFitter_h_INCLUDED
