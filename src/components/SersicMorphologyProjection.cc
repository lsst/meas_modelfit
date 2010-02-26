#include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
#include <ndarray/eigen.hpp>

namespace multifit = lsst::meas::multifit;
namespace components = multifit::components;

lsst::afw::geom::Extent2I 
components::SersicMorphologyProjection::getDimensions() const {
    return getKernelDimensions() * 4;
}

multifit::EllipseGridTransform::ConstPtr 
components::SersicMorphologyProjection::computeEllipseGridTransform() const {
    lsst::afw::geom::ellipses::BaseCore::Ptr transformedEllipse(
        getMorphology()->computeBoundingEllipseCore()->transform(*getTransform()).copy()
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

    //only one linear parameter: flux
    //grab the 0th element of the _linearParameterDerivative
    PartialDerivative output(_linearParameterDerivative[0]);

    InterpolationFunction::ConstPtr indexFunctor(
        getMorphology()->getSersicIndexFunctor()
    );
    lsst::afw::geom::LinearTransform egt = *computeEllipseGridTransform();    

    lsst::afw::geom::Point2D point(0);
    //outer loop over rows
    for (RowIter i(output.begin()), end(output.end()); i != end; ++i) {
        point.setX(0.0);
        //inner loop over pixels in each row
        for (PixIter j(i->begin()), rowEnd(i->end()); j != rowEnd; ++j) {
            double radius = egt(point).asVector().norm();
            *j = (*indexFunctor)(radius);
            point.setX(point.getX() + 1.0);
        }
        point.setY(point.getY() + 1.0);
    }

    return _linearParameterDerivative;
}

ndarray::FourierArray<multifit::Pixel,3,3> 
components::SersicMorphologyProjection::computeProjectedParameterDerivative() {
    typedef ndarray::Array<std::complex<Pixel>, 3> TransposedArray;
    typedef TransposedArray::Iterator RowIter;
    typedef TransposedArray::Reference::Iterator PixIter;

    //reinterpret the derivative array as (y, x, params)
    //allows for iterating over y, x
    TransposedArray output(
        _projectedParameterDerivative.getBase().transpose(ndarray::makeVector(1,2,0))
    );

    SersicMorphology::ConstPtr morphology(getMorphology());
    InterpolationFunction::ConstPtr indexFunctor(
        morphology->getSersicIndexFunctor()
    );
    EllipseGridTransform::ConstPtr ellipseGridPtr(computeEllipseGridTransform());
    EllipseGridTransform::DerivativeMatrix dEllipse(ellipseGridPtr->dEllipse());
    lsst::afw::geom::LinearTransform egt = *ellipseGridPtr;

    lsst::afw::geom::Point2D point(0);

    //outer loop over rows
    for (RowIter i(output.begin()), end(output.end()); i != end; ++i) {
        point.setX(0.0);
        //inner loop over pixels in each row
        for (PixIter j(i->begin()), rowEnd(i->end()); j != rowEnd; ++j) {
            //transform the point onto ellipse grid.
            lsst::afw::geom::Point2D ellipsePoint(egt(point));
            double radius = egt(ellipsePoint).asVector().norm();            
            //Use the row-functor over the sersic cache to compute the
            //partial derivative of the model in fourier space w.r.t to the
            //radius, then multiply by the partial derivative of the radius
            //w.r.t to the ellipse parameters (first 3 nonlinear model 
            //paramters)
            ndarray::viewAsEigen(*j).start<3>() = ( 
                indexFunctor->dParams(radius) * 
                (ellipsePoint.asVector()/radius).transpose() * 
                egt.dTransform(point) * dEllipse
            ).cast< std::complex<Pixel> >();
            
            //the last nonlinear parameter is the sersic index
            //Use the col-functor over the sersic cache to compute this last
            //partial derivative w.r.t sersic index
            InterpolationFunction::ConstPtr radiusFunctor(
                morphology->getCache()->getColFunctor(radius)
            );
            ndarray::viewAsEigen(*j).end<1>() << static_cast<std::complex<Pixel> >(
                radiusFunctor->dParams(morphology->getSersicIndex())
            );    
        }
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
    
    m->block<3,3>(0,0) << getMorphology()->computeBoundingEllipseCore()->transform(*getTransform()).d();
    (*m)(3,3) = 1;
    return ParameterJacobianMatrixPtr(m);
}

components::MorphologyProjection::TransformJacobianMatrixPtr
components::SersicMorphologyProjection::computeTransformParameterJacobian() const {
    TransformJacobianMatrix * m = new TransformJacobianMatrix(
        TransformJacobianMatrix::Zero(getNonlinearParameterSize())
    );    
    m->block<3,4>(0,0) << getMorphology()->computeBoundingEllipseCore()->transform(
        *getTransform()
    ).dTransform();
    return TransformJacobianMatrixPtr(m);
    
}

