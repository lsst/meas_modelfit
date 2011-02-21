// -*- LSST-C++ -*-
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

#ifndef LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletConvolution
#define LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletConvolution

#include "lsst/ndarray.h"
#include "lsst/afw/math/shapelets.h"
#include "lsst/afw/geom/ellipses.h"

#include <boost/scoped_ptr.hpp>

namespace lsst { namespace meas { namespace multifit { namespace shapelets {

/**
 *  @brief A parametrized matrix that performs a convolution in shapelet space.
 */
class ShapeletConvolution : private boost::noncopyable {
public:

    typedef boost::shared_ptr<ShapeletConvolution> Ptr;

    /**
     *  @brief Evaluate a shapelet convolution matrix in the given array.
     *
     *  @param[out]    matrix    A matrix that maps unconvolved shapelet coefficients (columns) to
     *                           convolved shapelet coefficients (rows).  Must already be allocated
     *                           to the correct shape.
     *  @param[in/out] ellipse   On input, the ellipse core of the unconvolved shapelet expansion.
     *                           On output, the ellipse core of the convolved shapelet expansion.
     */
    void evaluate(
        lsst::ndarray::Array<double,2,2> const & matrix, 
        lsst::afw::geom::ellipses::BaseCore & ellipse
    ) const;

    int getColOrder() const;

    int getRowOrder() const;

    lsst::ndarray::Vector<int,2> getShape() const {
        return ndarray::makeVector(
            afw::math::shapelets::computeSize(getRowOrder()),
            afw::math::shapelets::computeSize(getColOrder())
        );
    }

    ShapeletConvolution(
        int colOrder, 
        lsst::afw::math::shapelets::EllipticalShapeletFunction const & psf
    );

    ~ShapeletConvolution();

private:
    class Impl;

    boost::scoped_ptr<Impl> _impl;
};


}}}} // namespace lsst::meas::multifit::shapelets

#endif // !LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletConvolution
