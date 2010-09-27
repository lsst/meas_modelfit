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
 * Declaration of class PointSourceMorphologyProjection
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H

#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class PointSourceMorphology;

/**
 * Subclass of FourierMorphologyProjection used to represent projections of 
 * static point-sources
 */
class PointSourceMorphologyProjection : public FourierMorphologyProjection {
public:
    typedef boost::shared_ptr<PointSourceMorphologyProjection> Ptr;
    typedef boost::shared_ptr<PointSourceMorphologyProjection const> ConstPtr;
    typedef MorphologyProjection::ParameterJacobianMatrix ParameterJacobianMatrix;
    typedef MorphologyProjection::TransformJacobianMatrix TransformJacobianMatrix;
    typedef MorphologyProjection::ParameterJacobianMatrixPtr ParameterJacobianMatrixPtr;
    typedef MorphologyProjection::TransformJacobianMatrixPtr TransformJacobianMatrixPtr;
    
    PointSourceMorphologyProjection(
        boost::shared_ptr<PointSourceMorphology const> const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions, 
        lsst::afw::geom::AffineTransform const & transform
    ); 

    /// Imutable access to the PointSourceMorphology this is a projection of.
    boost::shared_ptr<PointSourceMorphology const> getMorphology() const { 
        return boost::static_pointer_cast<PointSourceMorphology const>(
            MorphologyProjection::getMorphology()
        );
    }

    // MorphologyProjection --------------------------------------------------
    virtual ParameterJacobianMatrixPtr computeProjectedParameterJacobian() const;
    virtual TransformJacobianMatrixPtr computeTransformParameterJacobian() const;

    // FourierMorphologyProjection --------------------------------------------
    virtual lsst::afw::geom::Extent2I getDimensions() const {
        return getKernelDimensions() + getPadding() * 2;
    }

    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative();
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative();
    
private:
    friend class PointSourceMorphology;

    ndarray::FourierArray<Pixel,3,3> _linearParameterDerivative;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H
