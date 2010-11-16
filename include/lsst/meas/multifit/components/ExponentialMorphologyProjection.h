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
 
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_EXPONENTIAL_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_EXPONENTIAL_MORPHOLOGY_PROJECTION_H

#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"
#include "lsst/meas/multifit/components/ExponentialMorphology.h"
#include "lsst/meas/multifit/EllipseGridTransform.h"
#include <ndarray/fft.hpp>
namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class ExponentialMorphologyProjection : public FourierMorphologyProjection {
public:
    typedef boost::shared_ptr<ExponentialMorphologyProjection> Ptr;
    typedef boost::shared_ptr<ExponentialMorphologyProjection const> ConstPtr;
    
    virtual lsst::afw::geom::Extent2I getDimensions() const;

    virtual ParameterJacobianMatrixPtr computeProjectedParameterJacobian() const;
    virtual TransformJacobianMatrixPtr computeTransformParameterJacobian() const; 

    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative();
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative();

    ExponentialMorphology::ConstPtr getMorphology() const {
        return boost::static_pointer_cast<ExponentialMorphology const>(
            MorphologyProjection::getMorphology()
        );
    }
protected:
    friend class ExponentialMorphology;
    typedef ndarray::FourierTransform<Pixel, 2> FFT;

    virtual void _handleLinearParameterChange();
    virtual void _handleNonlinearParameterChange();

    ExponentialMorphologyProjection(
        ExponentialMorphology::ConstPtr const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions,
        lsst::afw::geom::AffineTransform const & transform
    ) : FourierMorphologyProjection(morphology, kernelDimensions, transform), 
        _dimensions(),
        _validProducts(0)
    {
        _recomputeDimensions();
    }

    EllipseGridTransform::ConstPtr computeEllipseGridTransform() const;
    ndarray::FourierArray<Pixel, 3, 3> _projectedParameterDerivative;
    ndarray::FourierArray<Pixel, 3, 3> _linearParameterDerivative;

private:
    static double const UNCONVOLVED_RADIUS_FACTOR = 6;

    void _recomputeDimensions();
    lsst::afw::geom::Extent2I _dimensions;

    int _validProducts;
    enum Products {
        LINEAR_PARAMETER_DERIVATIVE=1,
        PROJECTED_PARAMETER_DERIVATIVE=2,
    };
};

}}}}
#endif
