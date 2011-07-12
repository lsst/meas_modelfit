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

    /// @brief Transform the matrix.
    void transform(lsst::afw::geom::AffineTransform const & t);

    /// @brief Transform the matrix with a simple scaling transform.
    void scale(double factor);

    /**
     *  @brief Modify an ellipse to include multipole moments that correspond to the given coefficients.
     */
    void applyMoments(
        lsst::afw::geom::ellipses::Ellipse & ellipse,
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients
    ) const;

    /// @brief Return the multipole matrix as an array.
    lsst::ndarray::Array<Pixel const,2,2> getArray() const { return _array; }

    /// @brief Construct from an array (always deep).
    explicit MultipoleMatrix(lsst::ndarray::Array<Pixel const,2,2> const & array) :
        _array(ndarray::copy(array))
    {}

    /// @brief Construct from an array (deep by default).
    explicit MultipoleMatrix(lsst::ndarray::Array<Pixel,2,2> const & array, bool deep=true) :
        _array(deep ? ndarray::copy(array) : array)
    {}

    /// @brief Copy construct (shallow by default).
    MultipoleMatrix(MultipoleMatrix const & other, bool deep=false) :
        _array(deep ? ndarray::copy(other._array) : other._array)
    {}

private:
    ndarray::Array<Pixel,2,2> _array;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MultipoleMatrix
