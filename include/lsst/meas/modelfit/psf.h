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

#ifndef LSST_MEAS_MODELFIT_psf_h_INCLUDED
#define LSST_MEAS_MODELFIT_psf_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/shapelet/FunctorKeys.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/modelfit/common.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  Control object used to define one piece of multishapelet fit to a PSF model; see PsfFitterControl
 */
class PsfFitterComponentControl {
public:

    PsfFitterComponentControl(int order_=0, double radiusFactor_=1.0, double radiusMin_=0.25) :
        order(order_), radiusFactor(radiusFactor_), radiusMin(radiusMin_)
    {}

    LSST_CONTROL_FIELD(
        order, int,
        "shapelet order for this component; negative to disable this component completely"
    );

    LSST_CONTROL_FIELD(
        radiusFactor, double,
        "Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either "
        "the second-moments radius of the PSF image (in an initial fit), or the radius of the primary "
        "component in a previous fit.  Ignored if the previous fit included this component (as then we "
        "can just use that radius)."
    );

    LSST_CONTROL_FIELD(
        radiusMin, double,
        "Require ellipses to be at least this large in all directions (in pixels)."
    );

};

/**
 * Control object used to configure a multishapelet fit to a PSF model; see PsfFitter.
 *
 *  The default configuration corresponds to fitting an elliptical double-Gaussian, in which each component
 *  can have different radii, positions, and ellipticities.  While the fitter can support much more complex
 *  models, at present, fitting these is prohibitively slow, and is not recommended in production
 *  environments.
 */
class PsfFitterControl {
public:

    PsfFitterControl() :
        inner(-1, 0.5), primary(0, 1.0), wings(0, 2.0), outer(-1, 4.0),
        maxIterations(20),
        absTol(1E-6),
        relTol(1E-6)
    {}

    LSST_NESTED_CONTROL_FIELD(
        inner, lsst.meas.modelfit.modelfitLib, PsfFitterComponentControl,
        "Innermost shapelet expansion, used to fit PSFs with very sharp cores"
    );

    LSST_NESTED_CONTROL_FIELD(
        primary, lsst.meas.modelfit.modelfitLib, PsfFitterComponentControl,
        "Primary shapelet expansion, typically used to fit the bulk of the PSF "
    );

    LSST_NESTED_CONTROL_FIELD(
        wings, lsst.meas.modelfit.modelfitLib, PsfFitterComponentControl,
        "Wing shapelet expansion (between primary and outer), typically used to fit the wings of the PSF"
    );

    LSST_NESTED_CONTROL_FIELD(
        outer, lsst.meas.modelfit.modelfitLib, PsfFitterComponentControl,
        "Outermost shapelet expansion, used to fit PSFs with very broad wings"
    );

    LSST_CONTROL_FIELD(
        maxIterations, int,
        "Maximum number of expectation-maximization iterations"
    );

    LSST_CONTROL_FIELD(
        absTol, double,
        "Absolute tolerance for expectation-maximization convergence: best-fit linear coefficients "
        "must not change more than this to declare convergence."
    );

    LSST_CONTROL_FIELD(
        relTol, double,
        "Relative tolerance for expectation-maximization convergence: best-fit linear coefficients "
        "must not change more than this to declare convergence (relative to mean of absolute values)."
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
 */
class PsfFitter {
public:

    /// Initialize the fitter class with the given control object.
    explicit PsfFitter(PsfFitterControl const & ctrl);

    /**
     *  Add fields to a Schema that can be used to store the MultiShapeletFunction returned by apply().
     *
     *  @param[in,out]   schema    Schema to add fields to.
     *  @param[in]       prefix    Field name prefix for all fields.
     *  @return a FunctorKey that can get/set MultiShapeletFunctions that match the configuration of this
     *          fitter on a record.
     */
    shapelet::MultiShapeletFunctionKey addModelFields(
        afw::table::Schema & schema,
        std::string const & prefix
    ) const;

    /**
     *  Adapt a differently-configured previous fit to be used as an starting point for this PsfFitter.
     *
     *  @param[in] previous          The return value of apply() from a differently-configured
     *                               instance of PsfFitter.
     *  @param[in] previousFitter    The PsfFitter used to create previous.
     *
     *  @return a new MultiShapelet function that may be passed directly to apply().  When possible,
     *  the ellipse and shapelet coefficeints will be copied from previousFit; higher-order coefficients
     *  will be set to zero, and any components used in this but unused in the previous fit will have their
     *  ellipses set relative to the previous fit's "primary" component.
     */
    shapelet::MultiShapeletFunction adapt(
        shapelet::MultiShapeletFunction const & previous,
        PsfFitter const & previousFitter
    ) const;

    /**
     *  Guess an initial model for the PSF shapelet approximation from just the moments of the PSF image.
    *
     *  @param[in]  moments     Second moments of the PSF, typically result of Psf::computeShape() or running
     *                          some other adaptive moments code on the PSF image.
     */
    shapelet::MultiShapeletFunction makeInitial(afw::geom::ellipses::Quadrupole const & moments) const;

    //@{
    /**
     *  Perform a fit to a PSF image, using a previous fit as a starting point
     *
     *  @param[in,out]  model   The model to fit, with configuration (number of components and orders)
     *                          already consistent with this fitter.  Also used to set initial ellipses
     *                          for the fit.
     *  @param[in]  image       The image to fit to, typically the result of Psf::computeKernelImage().  The
     *                          image's xy0 should be set such that the center of the PSF is at (0,0).
     *
     *  @return The number of iterations of the expectation-maximization algorithm.  If equal to
     *          maxIterations, the fitter did not converge.
     */
    int apply(
        shapelet::MultiShapeletFunction & model,
        afw::image::Image<Pixel> const & image
    ) const;
    int apply(
        shapelet::MultiShapeletFunction & model,
        afw::image::Image<double> const & image
    ) const {
        return apply(model, afw::image::Image<float>(image, true));
    }
    //@}

    int getComponentCount() const { return _ctrls.size(); }

    int getMaxIterations() const { return _maxIterations; }

    ~PsfFitter();

private:

    class Impl;

    void adaptComponent(
        shapelet::MultiShapeletFunction & current,
        shapelet::MultiShapeletFunction const & previous,
        PsfFitter const & previousFitter,
        int currentIndex,
        int previousIndex
    ) const;

    int _maxIterations;
    double _absTol;
    double _relTol;
    int _innerIndex;
    int _primaryIndex;
    int _wingsIndex;
    int _outerIndex;
    std::vector<PsfFitterComponentControl> _ctrls;
    mutable PTR(Impl) _impl;  // holds state specific to a certain PSF size; will be reset if that changes.
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_psf_h_INCLUDED
