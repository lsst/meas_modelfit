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
 
#include <Eigen/Geometry>
#include "lsst/meas/multifit/EllipseGridTransform.h"

namespace multifit = lsst::meas::multifit;

Eigen::Matrix2d const & multifit::EllipseGridTransform::_dA() {
    static Eigen::Matrix2d m = 
        (Eigen::Matrix2d() << 1.0, 0.0, 0.0, 0.0).finished();
    return m;
}

Eigen::Matrix2d const & multifit::EllipseGridTransform::_dB() {
    static Eigen::Matrix2d m = 
        (Eigen::Matrix2d() << 0.0, 0.0, 0.0, 1.0).finished();
    return m;
}

multifit::EllipseGridTransform::EllipseGridTransform(
    lsst::afw::geom::ellipses::BaseCore const & ellipse,
    lsst::afw::geom::ExtentI const & realDimensions
) : _axes(), _matrices() {
    boost::shared_ptr<Matrices> matrices(new Matrices());
    matrices->jacobian = _axes.dAssign(ellipse);
    matrices->rotation = Eigen::Rotation2D<double>(
        -_axes[lsst::afw::geom::ellipses::Axes::THETA]
    );
    matrices->factor << 
        (-2.0*M_PI / realDimensions.getX()), 0.0,
        0.0, (-2.0*M_PI / realDimensions.getY());
    matrices->scaling <<
        _axes[lsst::afw::geom::ellipses::Axes::A], 0.0,
        0.0, _axes[lsst::afw::geom::ellipses::Axes::B];
    matrices->tail = matrices->rotation * matrices->factor;
    _matrices = matrices;
}

multifit::EllipseGridTransform::DerivativeMatrix 
multifit::EllipseGridTransform::dEllipse() const {
    DerivativeMatrix m = DerivativeMatrix::Zero(
        DerivativeMatrix::RowsAtCompileTime, 
        DerivativeMatrix::ColsAtCompileTime
    );
    Eigen::Matrix2d dTheta(Eigen::Rotation2D<double>(
            -_axes[lsst::afw::geom::ellipses::Axes::THETA] - M_PI_2
        )
    );
    Eigen::Map<Eigen::Matrix2d> dT_dA(m.col(0).data());
    Eigen::Map<Eigen::Matrix2d> dT_dB(m.col(1).data());
    Eigen::Map<Eigen::Matrix2d> dT_dTheta(m.col(2).data());
    dT_dA = _dA() * _matrices->tail;
    dT_dB = _dB() * _matrices->tail;
    dT_dTheta = _matrices->scaling * dTheta * _matrices->factor;
    return m * _matrices->jacobian;
}
