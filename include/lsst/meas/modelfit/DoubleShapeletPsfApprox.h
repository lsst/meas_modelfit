// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 LSST/AURA.
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

#ifndef LSST_MEAS_MODELFIT_DoubleShapeletPsfApprox_h_INCLUDED
#define LSST_MEAS_MODELFIT_DoubleShapeletPsfApprox_h_INCLUDED

#include "lsst/afw/geom.h"
#include "lsst/shapelet/FunctorKeys.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/modelfit/optimizer.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  Control object used to configure a 2-shapelet fit to a PSF model; see DoubleShapeletPsfApproxAlgorithm.
 */
class DoubleShapeletPsfApproxControl {
public:

    DoubleShapeletPsfApproxControl() :
        innerOrder(2), outerOrder(1),
        radiusRatio(2.0), peakRatio(0.1),
        minRadius(1.0), minRadiusDiff(0.5), maxRadiusBoxFraction(0.4)
    {}

    LSST_CONTROL_FIELD(innerOrder, int, "Shapelet order of inner expansion (0 == Gaussian)");

    LSST_CONTROL_FIELD(outerOrder, int, "Shapelet order of outer expansion (0 == Gaussian)");

    LSST_CONTROL_FIELD(radiusRatio, double, "Initial outer radius divided by inner radius");

    LSST_CONTROL_FIELD(
        peakRatio, double,
        "Initial outer Gaussian peak height divided by inner Gaussian peak height"
    );

    LSST_CONTROL_FIELD(
        minRadius, double,
        "Don't allow the semi-minor radius of any component to drop below this value (pixels)"
    );

    LSST_CONTROL_FIELD(
        minRadiusDiff, double,
        "Don't allow the determinant radii of the two components to differ by less than this (pixels)"
    );

    LSST_CONTROL_FIELD(
        maxRadiusBoxFraction, double,
        "Don't allow the semi-major radius of any component to go above this fraction of the PSF image width"
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.modelfit.optimizer, OptimizerControl,
        "Configuration of the optimizer used by DoubleShapeletPsfsApproxAlgorithm::fitProfile()."
    );

};


/**
 *  An algorithm that fits a 2-component shapelet approximation to the PSF model.
 *
 *  This algorithm fits a similar model to the one that has been used on the
 *  HSC side for several data releases, but uses the improved optimizer in
 *  meas_modelfit instead of the one in (the now defunct)
 *  meas_extensions_multiShapelet and using moments to reduce the number of
 *  parameters in the fit.  It is faster and more robust than
 *  GeneralShapeletPsfApprox, but much less flexible.
 *
 *  The model is not a fully general one, even for two shapelet components. We
 *  tie the ellipticities of the two components together, hold the center
 *  fixed, and keep their radii fixed at a value set in the configuration;
 *  this means we fit only one set of ellipse parameters, instead of two.
 */
class DoubleShapeletPsfApproxAlgorithm : public meas::base::SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    static base::FlagDefinitionList const & getFlagDefinitions();
    static base::FlagDefinition const FAILURE;
    static base::FlagDefinition const INVALID_POINT_FOR_PSF;
    static base::FlagDefinition const INVALID_MOMENTS;
    static base::FlagDefinition const MAX_ITERATIONS;

    typedef DoubleShapeletPsfApproxControl Control;

    /// Failure modes passed by MeasurementErrors thrown by this class.

    /**
     *  Construct an algorithm for use with record outputs.
     *
     *  Essentially all of the functionality of this class is accessible through static methods,
     *  so an instance only needs to be constructed if you want to store the results in a record object.
     *
     *  @param[in]      ctrl      Control object specifying the details of the
     *                            model and how to fit it.
     *  @param[in]      name      Name of the algorithm; will be used as the
     *                            prefix for all field names.
     *  @param[in,out]  schema    Schema to which fields will be added.  Must
     *                            already contain a slot_Centroid alias
     *                            pointing to fields with position information.
     */
    DoubleShapeletPsfApproxAlgorithm(
        Control const & ctrl,
        std::string const & name,
        afw::table::Schema & schema
    );

    /**
     *  Create a MultiShapeletFunction with orders and radii and amplitude ratios from the control object.
     *
     *  This creates a circular two-component MultiShapeletFunction with unit
     *  total flux and unit circle moments, with zeroth-order amplitudes
     *  matching ctrl.peakRatio and ellipse radii matching ctrl.radiusRatio.
     */
    static shapelet::MultiShapeletFunction initializeResult(Control const & ctrl);

    /**
     *  Update a MultiShapeletFunction's ellipses to match the first and second moments of a PSF image.
     *
     *  The relative ampltidues, radii, and ellipticities of the ellipses
     *  will not be modified; both ellipses will simply be transformed
     *  (together) such that the moments of the MultiShapeletFunction match
     *  the (unweighted) moments of the image, and the two components will be
     *  scaled together to match the total flux in the image.
     *
     *  @param[in,out] result    MultiShapeletFunction to update.  Must have
     *                           unit circle moments.  Probably the result of
     *                           a call to initializeResult().
     *  @param[in]     ctrl      Control object specifying the details of the
     *                           model and how to fit it.
     *  @param[in]     psfImage  Image of the PSF to fit.  Should have xy0 set
     *                           such that (0,0) is at the center of the image.
     *                           Probably the result of a call to
     *                           Psf::computeKernelImage().
     *
     *  @throw meas::base::MeasurementError if the resulting moments are
     *         invalid, either because they have negative determinant or the
     *         radii are out of bounds.
     */
    static void fitMoments(
        shapelet::MultiShapeletFunction & result,
        Control const & ctrl,
        afw::image::Image<Scalar> const & psfImage
    );

    /**
     *  Return an Objective object that can be used to fit the profile of the model.
     *
     *  By "profile" we mean the two radii and two zeroth-order amplitudes,
     *  which we fit in fitProfile() while holding all other parameters fixed.
     *
     *  This function is probably not useful to most users, who should just
     *  call fitProfile instead; it's public mostly to make it easier to test
     *  functionality internal to fitProfile().
     *
     *  The parameters expected by the Objective are:
     *   - The zeroth-order amplitude of the inner component.
     *   - The zeroth-order amplitude of the outer component.
     *   - The determinant radius of the inner component relative to the determinant
     *     radius of the moments ellipse.
     *   - The determinant radius of the outer component relative to the determinant
     *     radius of the moments ellipse.
     *
     *  @param[in]   moments    Ellipse computed from the moments of the
     *                          MultiShapeletFunction to be fit.
     *  @param[in]   ctrl       Control object specifying the details of the
     *                          model and how to fit it.
     *  @param[in]     psfImage  Image of the PSF to fit.  Should have xy0 set
     *                           such that (0,0) is at the center of the image.
     *                           Probably the result of a call to
     *                           Psf::computeKernelImage().
     */
    static std::shared_ptr<OptimizerObjective> makeObjective(
        afw::geom::ellipses::Ellipse const & moments,
        Control const & ctrl,
        afw::image::Image<Scalar> const & psfImage
    );


    /**
     *  Update a MultiShapeletFunction's zeroth-order profile by fitting radii and amplitudes.
     *
     *  Holding the center positions and ellipticities fixed, the radii and amplitudes of the
     *  two Gaussians are fit simultaneously with an optimizer.  Higher-order terms will be
     *  ignored.
     *
     *  @param[in,out] result    MultiShapeletFunction to update.  Probably the
     *                           result of a call to fitMoments().
     *  @param[in]     ctrl      Control object specifying the details of the
     *                           model and how to fit it.
     *  @param[in]     psfImage  Image of the PSF to fit.  Should have xy0 set
     *                           such that (0,0) is at the center of the image.
     *                           Probably the result of a call to
     *                           Psf::computeKernelImage().
     */
    static void fitProfile(
        shapelet::MultiShapeletFunction & result,
        Control const & ctrl,
        afw::image::Image<Scalar> const & psfImage
    );

    /**
     *  Update a MultiShapeletFunction's higher-order shapelet terms, holding everything else fixed.
     *
     *  All ellipse parameters and the zeroth-order (Gaussian) shapelet coefficients will be held fixed.
     *
     *  @param[in,out] result    MultiShapeletFunction to update.  Probably the
     *                           result of a call to fitProfile().
     *  @param[in]     ctrl      Control object specifying the details of the
     *                           model and how to fit it.
     *  @param[in]     psfImage  Image of the PSF to fit.  Should have xy0 set
     *                           such that (0,0) is at the center of the image.
     *                           Probably the result of a call to
     *                           Psf::computeKernelImage().
     */
    static void fitShapelets(
        shapelet::MultiShapeletFunction & result,
        Control const & ctrl,
        afw::image::Image<Scalar> const & psfImage
    );

    /**
     *  Run all fitting stages on the Psf attached to the given Exposure, saving the results in measRecord.
     *
     *  We first call fitMoments(), then fitProfile(), then fitShapelets().
     */
    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    /**
     *  Handle failures caught by the measurement plugin system, setting failure flags as appropriate.
     */
    void fail(
        afw::table::SourceRecord & measRecord,
        lsst::meas::base::MeasurementError * error=nullptr
    ) const;

private:
    Control _ctrl;
    meas::base::SafeCentroidExtractor _centroidExtractor;
    shapelet::MultiShapeletFunctionKey _key;
    lsst::meas::base::FlagHandler _flagHandler;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_DoubleShapeletPsfApprox_h_INCLUDED
