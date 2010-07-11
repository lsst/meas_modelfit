// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
/**
 * @file
 * Declaration of class FourierModelProjection
 */
#ifndef LSST_MEAS_MULTIFIT_FOURIER_MODEL_PROEJECTION_H
#define LSST_MEAS_MULTIFIT_FOURIER_MODEL_PROEJECTION_H

#include <ndarray/fft_fwd.hpp>

#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/math/LocalKernel.h"
#include "lsst/afw/image/Utils.h"

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/WindowedFootprint.h"
#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 * A derived type of ComponentModelProjection which is convolved in fourier
 * space.
 */
class FourierModelProjection : public ComponentModelProjection {
public:
    typedef boost::shared_ptr<FourierModelProjection> Ptr;
    typedef boost::shared_ptr<FourierModelProjection const> ConstPtr;

    virtual int const getPsfParameterSize() const;

    virtual ~FourierModelProjection();

    /**
     * Immutable reference to the FourierMorphologyProjection this is based on
     */
    components::FourierMorphologyProjection::ConstPtr getMorphologyProjection() const {
        return boost::static_pointer_cast<components::FourierMorphologyProjection const>(
            ComponentModelProjection::getMorphologyProjection()
        );
    }

protected:
    friend class lsst::meas::multifit::ComponentModel;

    virtual void _convolve(PsfConstPtr const & psf);
    /// Determine if a valid PSF been provided 
    virtual bool isConvolved() const { return _localKernel; }

    virtual void _computeLinearParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);
    virtual void _computePsfParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);
    virtual void _computeTranslationDerivative(ndarray::Array<Pixel,2,1> const & matrix);
    virtual void _computeProjectedParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);

    virtual bool hasPsfParameterDerivative() const {
        return _localKernel->hasDerivatives() && _localKernel->getNParameters() > 0;
    }

    virtual void _handleLinearParameterChange();
    virtual void _handleNonlinearParameterChange();

    /**
     * Mutable referece to the FourierMorphologyProjection this is based on
     */
    components::FourierMorphologyProjection::Ptr _getMorphologyProjection() {
        return boost::static_pointer_cast<components::FourierMorphologyProjection>(
            ComponentModelProjection::_getMorphologyProjection()
        );
    }
private:
    typedef ndarray::FourierTransform<Pixel,2> FFT;
    
    FourierModelProjection(
        ComponentModel::ConstPtr const & model,
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint
    );

    void _setDimensions();

    void _applyKernel(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
        ndarray::FourierArray<Pixel const,2,2> const & kernel
    ) const;

    void _applyKernel(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end
    ) const;

    // PSF/Kernel information.
    lsst::afw::math::FourierLocalKernel::Ptr _localKernel; 
    // maps footprint to handler output arrays
    WindowedFootprint::Ptr _wf;     
    // bounding box of padded arrays relative to exposure
    lsst::afw::geom::BoxI _outerBBox; 
    // bounding box of unpadded arrays relative to _outerBBox
    lsst::afw::geom::BoxI _innerBBox;

    class Shifter;
    class LinearMatrixHandler;
    class NonlinearMatrixHandler;
    class PsfMatrixHandler;

    boost::scoped_ptr<Shifter> _shifter;
    boost::scoped_ptr<LinearMatrixHandler> _linearMatrixHandler;
    boost::scoped_ptr<NonlinearMatrixHandler> _nonlinearMatrixHandler;
    boost::scoped_ptr<PsfMatrixHandler> _psfMatrixHandler;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_FOURIER_MODEL_PROEJECTION_H

