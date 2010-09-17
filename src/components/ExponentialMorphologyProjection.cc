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
 
#include "lsst/meas/multifit/components/ExponentialMorphologyProjection.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/pex/exceptions/Runtime.h"
#include <ndarray/eigen.hpp>

namespace multifit = lsst::meas::multifit;
namespace components = multifit::components;

void components::ExponentialMorphologyProjection::_recomputeDimensions() {
    lsst::afw::geom::ellipses::BaseCore::Ptr transformedEllipse(
        getMorphology()->computeBoundingEllipseCore()->transform(*getTransform()).copy()
    );
    lsst::afw::geom::Extent2D ellipseBounds = transformedEllipse->computeDimensions();
    //grow the ellipse dimensions by a constant factor
    ellipseBounds *= UNCONVOLVED_RADIUS_FACTOR;

    lsst::afw::geom::Extent2I dimensions = lsst::afw::geom::makeExtentI(
        static_cast<int>(std::ceil(ellipseBounds.getX())),
        static_cast<int>(std::ceil(ellipseBounds.getY()))
    );

    dimensions += getKernelDimensions() + getPadding()*2;

    if(lsst::afw::geom::any(_dimensions.lt(dimensions))) {
        //if the current _dimensions are smaller than the newly computed
        //dimensions, then set _dimensions to be the newly computed dimensions
        //grown by some factor to reduce need for future resizing
        _dimensions = dimensions*1.5;

        //reallocate product arrays
        ndarray::shallow(_projectedParameterDerivative) = FFT::initializeK(
            ndarray::makeVector(
                getNonlinearParameterSize(), 
                _dimensions.getY(), 
                _dimensions.getX()
            )
        );
        ndarray::shallow(_linearParameterDerivative) = FFT::initializeK(
            ndarray::makeVector(
                getLinearParameterSize(), 
                _dimensions.getY(), 
                _dimensions.getX()
            )
        );   
    }
}

lsst::afw::geom::Extent2I 
components::ExponentialMorphologyProjection::getDimensions() const {
    return _dimensions;
}

multifit::EllipseGridTransform::ConstPtr 
components::ExponentialMorphologyProjection::computeEllipseGridTransform() const {
    lsst::afw::geom::ellipses::BaseCore::Ptr transformedEllipse(
        getMorphology()->computeBoundingEllipseCore()->transform(*getTransform()).copy()
    );
    return boost::make_shared<EllipseGridTransform>(
        *transformedEllipse, getDimensions()
    );
}
ndarray::FourierArray<multifit::Pixel,3,3> 
components::ExponentialMorphologyProjection::computeLinearParameterDerivative() {
    typedef ndarray::FourierArray<Pixel, 3, 3>::Reference PartialDerivative;
    typedef PartialDerivative::Iterator RowIter;
    typedef PartialDerivative::Reference::Iterator PixIter;

    if((_validProducts & LINEAR_PARAMETER_DERIVATIVE) == 0) {
        //only one linear parameter: flux
        //grab the 0th element of the _linearParameterDerivative
        PartialDerivative output(_linearParameterDerivative[0]);

        lsst::afw::geom::LinearTransform egt = *computeEllipseGridTransform();
        lsst::afw::geom::Extent2I dimensions = getDimensions();
        int midY = dimensions.getY()/2;
        lsst::afw::geom::Point2D point;
        int y=0, x;
        //outer loop over rows
        for (RowIter i(output.begin()), end(output.end()); i != end; ++i, ++y) {
            point.setY((y > midY) ? (y - dimensions.getY()) : y);
            x = 0;
            //inner loop over pixels in each row
            for (PixIter j(i->begin()), rowEnd(i->end()); j != rowEnd; ++j, ++x) {            
                point.setX(x);
                double k = egt(point).asVector().norm();
                *j = std::pow(1.0 + k*k, -1.5);
            }
        }
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }
    return _linearParameterDerivative;
}

ndarray::FourierArray<multifit::Pixel,3,3> 
components::ExponentialMorphologyProjection::computeProjectedParameterDerivative() {
    typedef ndarray::Array<std::complex<Pixel>, 3> TransposedArray;
    typedef TransposedArray::Iterator RowIter;
    typedef TransposedArray::Reference::Iterator PixIter;

    if((_validProducts & PROJECTED_PARAMETER_DERIVATIVE) == 0) {
        //reinterpret the derivative array as (y, x, params)
        //allows for iterating over y, x
        TransposedArray output(
            _projectedParameterDerivative.getBase().transpose(ndarray::makeVector(1,2,0))
        );

        ExponentialMorphology::ConstPtr morphology(getMorphology());
        EllipseGridTransform::ConstPtr ellipseGridPtr(computeEllipseGridTransform());
        EllipseGridTransform::DerivativeMatrix dEllipse(ellipseGridPtr->dEllipse());
        lsst::afw::geom::LinearTransform egt = *ellipseGridPtr;
        lsst::afw::geom::Extent2I dimensions = getDimensions();
        int midY = dimensions.getY()/2;
        lsst::afw::geom::Point2D point(0);

        int y=0,x;
        double dParams;
        //outer loop over rows
        for (RowIter i(output.begin()), end(output.end()); i != end; ++i, ++y) {
            point.setY((y>midY) ? (y-dimensions.getY()) : y); 
            x= 0;
            //inner loop over pixels in each row
            for (PixIter j(i->begin()), rowEnd(i->end()); j != rowEnd; ++j, ++x) {            
                point.setX(x);            
                //transform the point onto ellipse grid.
                lsst::afw::geom::Point2D ellipsePoint(egt(point));       
                double k = ellipsePoint.asVector().norm();  

                dParams = -3.0 * k * std::pow(1.0 + k*k, -2.5);
                //Use the row-functor over the exponential cache to compute the
                //partial derivative of the model in fourier space w.r.t to the
                //radius, then multiply by the partial derivative of the radius
                //w.r.t to the ellipse parameters (first 3 nonlinear model 
                //paramters)
                if(k < 1E-16) {
                    ellipsePoint.setX(1.0);
                    ellipsePoint.setY(1.0);
                } else {
                    ellipsePoint.setX(ellipsePoint.getX()/k);
                    ellipsePoint.setY(ellipsePoint.getY()/k);
                }
                
                ndarray::viewAsEigen(*j).start<3>() = ( 
                    dParams * ellipsePoint.asVector().transpose() * 
                    egt.dTransform(point) * dEllipse
                ).cast< std::complex<Pixel> >();
            }        
        }
        _validProducts |= PROJECTED_PARAMETER_DERIVATIVE;
    }

    return _projectedParameterDerivative;
}
    
components::MorphologyProjection::ParameterJacobianMatrixPtr
components::ExponentialMorphologyProjection::computeProjectedParameterJacobian() const { 
    ParameterJacobianMatrix * m = new ParameterJacobianMatrix(
        ParameterJacobianMatrix::Zero(
            getNonlinearParameterSize(),
            getNonlinearParameterSize()
        )
    );
    
    m->block<3,3>(0,0) << getMorphology()->computeBoundingEllipseCore()->transform(*getTransform()).d();
    return ParameterJacobianMatrixPtr(m);
}

components::MorphologyProjection::TransformJacobianMatrixPtr
components::ExponentialMorphologyProjection::computeTransformParameterJacobian() const {
    TransformJacobianMatrix * m = new TransformJacobianMatrix(
        TransformJacobianMatrix::Zero(getNonlinearParameterSize(), TransformJacobianMatrix::ColsAtCompileTime)
    );    
    m->block<3,4>(0,0) << getMorphology()->computeBoundingEllipseCore()->transform(
        *getTransform()
    ).dTransform();
    return TransformJacobianMatrixPtr(m);
    
}

void components::ExponentialMorphologyProjection::_handleNonlinearParameterChange() {
    _recomputeDimensions();

    _validProducts &= ~LINEAR_PARAMETER_DERIVATIVE;
    _validProducts &= ~PROJECTED_PARAMETER_DERIVATIVE;
}

void components::ExponentialMorphologyProjection::_handleLinearParameterChange() {
    _validProducts &= ~LINEAR_PARAMETER_DERIVATIVE;
}


