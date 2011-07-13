// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_MultipoleMatrix
#define LSST_MEAS_MULTIFIT_MultipoleMatrix

#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

class MultipoleMatrix {
public:

    enum { I0=0, IX, IY, IXX, IYY, IXY };

    lsst::afw::geom::Point2D computeCentroid(
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
        lsst::afw::geom::ellipses::Ellipse const & ellipse
    ) const;

    lsst::afw::geom::Point2D computeCentroid(
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients
    ) const;

    lsst::afw::geom::ellipses::Quadrupole computeQuadrupole(
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
        lsst::afw::geom::ellipses::Ellipse const & ellipse,
        lsst::afw::geom::Point2D const & centroid
    ) const;

    lsst::afw::geom::ellipses::Quadrupole computeQuadrupole(
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
        lsst::afw::geom::Point2D const & centroid
    ) const;

    lsst::afw::geom::ellipses::Ellipse computeEllipse(
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
        lsst::afw::geom::ellipses::Ellipse const & ellipse
    ) const;

    lsst::afw::geom::ellipses::Ellipse computeEllipse(
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients
    ) const;

    /// @brief Return the multipole matrix as an array.
    lsst::ndarray::Array<Pixel const,2,2> getArray() const { return _array; }

    /**
     *  @brief Update a grid parameter vector by modifying all centers and
     *         ellipses to include the moments defined by a coefficient vector.
     */
    static void updateParameters(
        PTR(grid::Grid) const & grid, 
        lsst::ndarray::Array<double const,1,1> const & parametersIn,
        lsst::ndarray::Array<double,1,1> const & parametersOut,
        lsst::ndarray::Array<Pixel const,1,1> & coefficients
    );

    /// @brief Construct from an array (shallow).
    explicit MultipoleMatrix(lsst::ndarray::Array<Pixel const,2,2> const & array) : _array(array) {}

    /// @brief Copy construct (shallow).
    MultipoleMatrix(MultipoleMatrix const & other) : _array(other._array) {}

private:
    ndarray::Array<Pixel const,2,2> _array;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MultipoleMatrix
