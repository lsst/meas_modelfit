// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_constants_h_INCLUDED
#define LSST_MEAS_MULTIFIT_constants_h_INCLUDED

#include "Eigen/Core"
#include "ndarray_fwd.h"

namespace lsst { namespace meas { namespace multifit {

//@{
/** Typedefs to be used for pixel and probability values throughout meas_multifit.
 *
 *  While it's perfectly sufficient to store pretty much everything in single-precision,
 *  there are some calculations we do in double-precision, particular sums over probability
 *  samples.  In those cases, we'll explicitly use double along with the usual Eigen
 *  matrix and vector typedefs.
 */
typedef float Pixel;
typedef Eigen::Matrix<Pixel,Eigen::Dynamic,Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<Pixel,Eigen::Dynamic,1> Vector;
typedef ndarray::Array<Pixel,1,1> Array1;
typedef ndarray::Array<Pixel,2,-1> Array2CM;
typedef ndarray::Array<Pixel,2,1> Array2RM;
//@}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_constants_h_INCLUDED
