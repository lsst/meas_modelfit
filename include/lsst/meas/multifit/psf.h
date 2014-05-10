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

#ifndef LSST_MEAS_MULTIFIT_psf_h_INCLUDED
#define LSST_MEAS_MULTIFIT_psf_h_INCLUDED

#include "boost/scoped_ptr.hpp"

#include "lsst/pex/config.h"

#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/Likelihood.h"

namespace lsst { namespace meas { namespace multifit {

class PsfFitterComponentControl {
public:

    PsfFitterComponentControl(int order_=0, double radiusFactor_=1.0) :
        order(order_), positionPriorSigma(0.0), ellipticityPriorSigma(0.0),
        radiusFactor(radiusFactor_), radiusPriorSigma(0.0)
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

class PsfFitterControl {
public:

    PsfFitterControl() : inner(-1, 0.5), primary(0, 1.0), wings(0, 2.0), outer(-1, 4.0) {}

    LSST_NESTED_CONTROL_FIELD(
        inner, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Innermost shapelet expansion, used to fit PSFs with very sharp cores"
    );

    LSST_NESTED_CONTROL_FIELD(
        primary, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Primary shapelet expansion, typically used to fit the bulk of the PSF "
    );

    LSST_NESTED_CONTROL_FIELD(
        wings, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Wing shapelet expansion (between primary and outer), typically used to fit the wings of the PSF"
    );

    LSST_NESTED_CONTROL_FIELD(
        outer, lsst.meas.multifit.multifitLib, PsfFitterComponentControl,
        "Outermost shapelet expansion, used to fit PSFs with very broad wings"
    );

};

class PsfFitter {
public:

    PsfFitter(PsfFitterControl const & ctrl);

    PTR(Model) getModel() const { return _model; }

    PTR(Prior) getPrior() const { return _prior; }

    shapelet::MultiShapeletFunction adapt(
        shapelet::MultiShapeletFunction const & previousFit,
        PTR(Model) previousModel
    ) const;

    shapelet::MultiShapeletFunction apply(
        afw::image::Image<Pixel> const & image,
        Scalar noiseSigma,
        afw::geom::ellipses::Quadrupole const & moments
    ) const;

    shapelet::MultiShapeletFunction apply(
        afw::image::Image<Pixel> const & image,
        Scalar noiseSigma,
        shapelet::MultiShapeletFunction const & initial
    ) const;

private:
    PsfFitterControl _ctrl;
    PTR(Model) _model;
    PTR(Prior) _prior;
};


class MultiShapeletPsfLikelihood : public Likelihood {
public:

    MultiShapeletPsfLikelihood(
        ndarray::Array<Pixel const,2,2> const & image,
        afw::geom::Point2I const & xy0,
        PTR(Model) model,
        Scalar sigma,
        ndarray::Array<Scalar const,1,1> const & fixed
    );

    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        bool doApplyWeights=true
    ) const;

    virtual ~MultiShapeletPsfLikelihood();

private:
    class Impl;
    boost::scoped_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_psf_h_INCLUDED
