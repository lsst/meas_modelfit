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
 
#include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
#include "lsst/meas/multifit/SersicCache.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/pex/exceptions/Runtime.h"
#include <ndarray/eigen.hpp>

namespace multifit = lsst::meas::multifit;
namespace components = multifit::components;

void components::SersicMorphologyProjection::_recomputeDimensions() {
    lsst::afw::geom::ellipses::BaseCore::Ptr transformedEllipse(
        getMorphology()->computeBoundingEllipseCore()->transform(getTransform()).copy()
    );
    lsst::afw::geom::Extent2D ellipseBounds = transformedEllipse->computeDimensions();
    //grow the ellipse dimensions by a constant factor
    ellipseBounds *= 3;

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
components::SersicMorphologyProjection::getDimensions() const {
    return _dimensions;
}

multifit::EllipseGridTransform::ConstPtr 
components::SersicMorphologyProjection::computeEllipseGridTransform() const {
    lsst::afw::geom::ellipses::BaseCore::Ptr transformedEllipse(
        getMorphology()->computeBoundingEllipseCore()->transform(getTransform()).copy()
    );
    return boost::make_shared<EllipseGridTransform>(
        *transformedEllipse, getDimensions()
    );
}
ndarray::FourierArray<multifit::Pixel,3,3> 
components::SersicMorphologyProjection::computeLinearParameterDerivative() {
    typedef ndarray::FourierArray<Pixel, 3, 3>::Reference PartialDerivative;
    typedef PartialDerivative::Iterator RowIter;
    typedef PartialDerivative::Reference::Iterator PixIter;

    if((_validProducts & LINEAR_PARAMETER_DERIVATIVE) == 0) {
        //only one linear parameter: flux
        //grab the 0th element of the _linearParameterDerivative
        PartialDerivative output(_linearParameterDerivative[0]);

        SersicCache::Interpolator::ConstPtr interpolator(getMorphology()->getInterpolator());
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
                try {
                    *j = (*interpolator)(k);
                }
                catch (lsst::pex::exceptions::InvalidParameterException & ) {
                    throw LSST_EXCEPT(
                        lsst::pex::exceptions::LogicErrorException,
                        (boost::format("k %1% value out of range. Check SersicCachePolicy")%k).str()
                    );
                }
            }
        }
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }
    return _linearParameterDerivative;
}

ndarray::FourierArray<multifit::Pixel,3,3> 
components::SersicMorphologyProjection::computeProjectedParameterDerivative() {
    typedef ndarray::Array<std::complex<Pixel>, 3> TransposedArray;
    typedef TransposedArray::Iterator RowIter;
    typedef TransposedArray::Reference::Iterator PixIter;

    if((_validProducts & PROJECTED_PARAMETER_DERIVATIVE) == 0) {
        //reinterpret the derivative array as (y, x, params)
        //allows for iterating over y, x
        TransposedArray output(
            _projectedParameterDerivative.getBase().transpose(ndarray::makeVector(1,2,0))
        );

        SersicMorphology::ConstPtr morphology(getMorphology());
        SersicCache::Interpolator::ConstPtr interpolator(morphology->getInterpolator());
        SersicCache::Interpolator::ConstPtr derivativeInterpolator(morphology->getDerivativeInterpolator());

        EllipseGridTransform::ConstPtr ellipseGridPtr(computeEllipseGridTransform());
        EllipseGridTransform::DerivativeMatrix dEllipse(ellipseGridPtr->dEllipse());
        lsst::afw::geom::LinearTransform egt = *ellipseGridPtr;
        lsst::afw::geom::Extent2I dimensions = getDimensions();
        int midY = dimensions.getY()/2;
        lsst::afw::geom::Point2D point(0);

        int y=0,x;
        double dK, dN;
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
                try {
                    dK = interpolator->d(k);
                } 
                catch (lsst::pex::exceptions::InvalidParameterException &) {
                    throw LSST_EXCEPT(
                        lsst::pex::exceptions::LogicErrorException,
                        (boost::format("k %1% value out of range. Check SersicCachePolicy")%k).str()
                    );

                }
                //Use the row-functor over the sersic cache to compute the
                //partial derivative of the model in fourier space w.r.t to the
                //radius, then multiply by the partial derivative of the radius
                //w.r.t to the ellipse parameters (first 3 nonlinear model 
                //paramters)
                if(k == 0.0) {
                    ellipsePoint.setX(1.0);
                    ellipsePoint.setY(1.0);
                } else {
                    ellipsePoint.setX(ellipsePoint.getX()/k);
                    ellipsePoint.setY(ellipsePoint.getY()/k);
                }
                
                ndarray::viewAsEigen(*j).start<3>() = ( 
                    dK * ellipsePoint.asVector().transpose() * 
                    egt.dTransform(point) * dEllipse
                ).cast< std::complex<Pixel> >();
                
                //the last nonlinear parameter is the sersic index
                //Use the col-functor over the sersic cache to compute this last
                //partial derivative w.r.t sersic index
                try {
                    dN = (*derivativeInterpolator)(k);
                } catch(lsst::pex::exceptions::InvalidParameterException &) {
                    throw LSST_EXCEPT(
                        lsst::pex::exceptions::LogicErrorException,
                        (boost::format("k %1% value out of range. Check SersicCachePolicy")%k).str()
                    );
                }
                ndarray::viewAsEigen(*j).end<1>() << static_cast<std::complex<Pixel> >(dN);    
            }        
        }
        _validProducts |= PROJECTED_PARAMETER_DERIVATIVE;
    }

    return _projectedParameterDerivative;
}
    
components::MorphologyProjection::ParameterJacobianMatrixPtr
components::SersicMorphologyProjection::computeProjectedParameterJacobian() const { 
    ParameterJacobianMatrix * m = new ParameterJacobianMatrix(
        ParameterJacobianMatrix::Zero(
            getNonlinearParameterSize(),
            getNonlinearParameterSize()
        )
    );
    
    m->block<3,3>(0,0) << getMorphology()->computeBoundingEllipseCore()->transform(getTransform()).d();
    (*m)(3,3) = 1;
    return ParameterJacobianMatrixPtr(m);
}

components::MorphologyProjection::TransformJacobianMatrixPtr
components::SersicMorphologyProjection::computeTransformParameterJacobian() const {
    TransformJacobianMatrix * m = new TransformJacobianMatrix(
        TransformJacobianMatrix::Zero(getNonlinearParameterSize(), TransformJacobianMatrix::ColsAtCompileTime)
    );    
    m->block<3,4>(0,0) << getMorphology()->computeBoundingEllipseCore()->transform(
        getTransform()
    ).dTransform();
    return TransformJacobianMatrixPtr(m);
    
}

void components::SersicMorphologyProjection::_handleNonlinearParameterChange() {
    _recomputeDimensions();

    _validProducts &= ~LINEAR_PARAMETER_DERIVATIVE;
    _validProducts &= ~PROJECTED_PARAMETER_DERIVATIVE;
}

void components::SersicMorphologyProjection::_handleLinearParameterChange() {
    _validProducts &= ~LINEAR_PARAMETER_DERIVATIVE;
}


