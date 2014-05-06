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

    PsfFitterComponentControl(int initialOrder_=0, int finalOrder_=0, double fiducialRadius_=1.0) :
        initialOrder(initialOrder_), finalOrder(finalOrder_),
        positionPriorSigma(0.0), ellipticityPriorSigma(0.0),
        fiducialRadius(fiducialRadius_), radiusPriorSigma(0.0)
    {}

    LSST_CONTROL_FIELD(
        initialOrder, int,
        "shapelet order in initial fit; -1 for none"
    );
    LSST_CONTROL_FIELD(
        finalOrder, int,
        "shapelet order in final fit; -1 for none"
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
        fiducialRadius, double,
        "fiducial radius of the shapelet expansion, in units of the 2nd moments radius of the image"
    );
    LSST_CONTROL_FIELD(
        radiusPriorSigma, double,
        "sigma in a Gaussian prior on ln(radius/fiducialRadius)"
    );

};

class PsfFitterControl {
public:

    PsfFitterControl() : inner(-1, 0, 0.5), primary(0, 5, 1.0), wings(0, 0, 2.0), outer(-1, 0, 4.0) {}

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

    PTR(Model) getInitialModel() const;

    PTR(Model) getFinalModel() const;

    PTR(Prior) getInitialPrior() const;

    PTR(Prior) getFinalPrior() const;

    shapelet::MultiShapeletFunction fitInitial(
        afw::image::Image<Pixel> const & image,
        Scalar noiseSigma
    );

    shapelet::MultiShapeletFunction fitFinal(
        afw::image::Image<Pixel> const & image,
        Scalar noiseSigma,
        shapelet::MultiShapeletFunction const & initialResult
    );

private:
    PsfFitterControl _ctrl;
    PTR(Model) _initialModel;
    PTR(Model) _finalModel;
    PTR(Prior) _initialPrior;
    PTR(Prior) _finalPrior;
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
