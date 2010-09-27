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
 * Implementation of class PointSourceMorphologyProjection
 */
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"

namespace multifit = lsst::meas::multifit;
namespace components = multifit::components;

components::PointSourceMorphologyProjection::ParameterJacobianMatrixPtr 
components::PointSourceMorphologyProjection::computeProjectedParameterJacobian() const {
    static const ParameterJacobianMatrixPtr m(new ParameterJacobianMatrix()); 
    // matrix has zero size
    return m;
}

components::PointSourceMorphologyProjection::TransformJacobianMatrixPtr
components::PointSourceMorphologyProjection::computeTransformParameterJacobian() const {
    static const TransformJacobianMatrixPtr m(new TransformJacobianMatrix());
    // matrix has zero size
    return m;
}

ndarray::FourierArray<multifit::Pixel,3,3>
components::PointSourceMorphologyProjection::computeLinearParameterDerivative() {
    return _linearParameterDerivative;

}

ndarray::FourierArray<multifit::Pixel,3,3>
components::PointSourceMorphologyProjection::computeProjectedParameterDerivative() {
    return ndarray::FourierArray<Pixel,3,3>();
}

/**
 * Construct a PointSourceMorphologyProjection.
 */
lsst::meas::multifit::components::PointSourceMorphologyProjection::PointSourceMorphologyProjection(
    PointSourceMorphology::ConstPtr const & morphology,
    lsst::afw::geom::Extent2I const & kernelDimensions, 
    lsst::afw::geom::AffineTransform const & transform
) : FourierMorphologyProjection(morphology,kernelDimensions,transform),
    _linearParameterDerivative()
{
    ndarray::shallow(_linearParameterDerivative) = ndarray::FourierTransform<Pixel,2>::initializeK(
        ndarray::makeVector(1, getDimensions().getY(), getDimensions().getX())
    );
    _linearParameterDerivative = 1.0;
}
