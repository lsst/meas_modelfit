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
#include "boost/make_shared.hpp"
#include "lsst/meas/multifit/EllipseGridTransform.h"

namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace ellipses = lsst::afw::geom::ellipses;

multifit::EllipseGridTransform::EllipseGridTransform(
    ellipses::BaseCore const & ellipse,
    geom::ExtentI const & realDimensions
) : _logShear(), _matrices(new Matrices) {
    // Do not use make_shared in the constructor for _matrices, as it bypasses
    // Eigen's special 128-bit aligned allocator.
    _matrices->jacobian = _logShear.dAssign(ellipse);
    Eigen::Matrix2d mI = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d mZ; mZ << 1.0, 0.0, 0.0,-1.0; // Z, X denote pauli spin matrices
    Eigen::Matrix2d mX; mX << 0.0, 1.0, 1.0, 0.0;
    double expk = std::exp(_logShear[ellipses::LogShear::KAPPA]);
    double gamma = _logShear.getGamma();
    double gamma1 = _logShear[ellipses::LogShear::GAMMA1];
    double gamma2 = _logShear[ellipses::LogShear::GAMMA2];
    double coshg = std::cosh(gamma);
    double sinhg = std::sinh(gamma);
    if (gamma > 1E-16) {
        double sinhg_g = sinhg / gamma;
        double f = (coshg - sinhg_g) / (gamma * gamma);
        _matrices->primary = expk * (mI * coshg + (gamma1 * mZ + gamma2 * mX) * sinhg_g);
        _matrices->dgamma1 = expk * (mI * gamma1 * sinhg_g + mZ * sinhg_g
                                     + (gamma1 * mZ + gamma2 * mX) * gamma1 * f);
        _matrices->dgamma2 = expk * (mI * gamma2 * sinhg_g + mX * sinhg_g
                                     + (gamma1 * mZ + gamma2 * mX) * gamma2 * f);
    } else {
        // TODO(?): carry these out to second order in the Taylor expansion
        _matrices->primary = expk * (mI * coshg + (gamma1 * mZ + gamma2 * mX));
        _matrices->dgamma1 = expk * mZ;
        _matrices->dgamma2 = expk * mX;
    }
    _matrices->tail << 
        (2.0*M_PI / realDimensions.getX()), 0.0,
        0.0, (2.0*M_PI / realDimensions.getY());
}

multifit::EllipseGridTransform::DerivativeMatrix 
multifit::EllipseGridTransform::dEllipse() const {
    DerivativeMatrix m = DerivativeMatrix::Zero();
    geom::LinearTransform dgamma1(_matrices->dgamma1 * _matrices->tail);
    geom::LinearTransform dgamma2(_matrices->dgamma1 * _matrices->tail);
    geom::LinearTransform dkappa(_matrices->primary * _matrices->tail);
    for (int n = 0; n < m.rows(); ++n) {
        m(n, 0) = dgamma1[n];
        m(n, 1) = dgamma2[n];
        m(n, 2) = dkappa[n];
    }
    return m * _matrices->jacobian;
}
