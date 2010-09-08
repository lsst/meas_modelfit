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
 
#ifndef LSST_MEAS_MULTIFIT_ELLIPSE_GRID_TRANSFORM_H
#define LSST_MEAS_MULTIFIT_ELLIPSE_GRID_TRANSFORM_H

#include <Eigen/Core>
#include "lsst/afw/geom.h"


namespace lsst {
namespace meas {
namespace multifit {

/**
 *  A differentiable expression object that computes the LinearTransform that, 
 *  when applied to a Fourier-space coordinate grid, maps the unit circle onto 
 *  an ellipse in real-space.
 */
class EllipseGridTransform {
public:
    typedef boost::shared_ptr<EllipseGridTransform> Ptr;
    typedef boost::shared_ptr<EllipseGridTransform const> ConstPtr;

    typedef Eigen::Matrix<double,4,3> DerivativeMatrix;
    
    explicit EllipseGridTransform(
        lsst::afw::geom::ellipses::BaseCore const & ellipse,
        lsst::afw::geom::ExtentI const & realDimensions
    );

    operator lsst::afw::geom::LinearTransform () const {
        return lsst::afw::geom::LinearTransform(_matrices->primary * _matrices->tail);
    }

    DerivativeMatrix dEllipse() const;

private:
    struct Matrices {
        Eigen::Matrix2d primary;
        Eigen::Matrix2d dgamma1;
        Eigen::Matrix2d dgamma2;
        Eigen::Matrix2d tail;
        Eigen::Matrix3d jacobian;
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW; 
    };

    lsst::afw::geom::ellipses::LogShear _logShear;
    boost::shared_ptr<Matrices> _matrices;

};

}}} // namespace lsst::meas::multifit

#endif 
