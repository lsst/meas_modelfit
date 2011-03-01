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

#ifndef LSST_MEAS_MULTIFIT_SAMPLING_ImportanceSample
#define LSST_MEAS_MULTIFIT_SAMPLING_ImportanceSample

#include "lsst/ndarray.h"
#include "lsst/ndarray/tables.h"

#include <boost/fusion/container/vector.hpp>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

struct ImportanceSample {
    typedef boost::fusion::vector< 
        ndarray::tables::Field<double>,
        ndarray::tables::Field<double>,
        ndarray::tables::Field<double>,
        ndarray::tables::Field<double,1>,
        ndarray::tables::Field<double,1>,
        ndarray::tables::Field<double,2>
        > FieldSequence;

    static ndarray::tables::Index<0> const PROPOSAL;
    static ndarray::tables::Index<1> const TARGET;
    static ndarray::tables::Index<2> const WEIGHT;
    static ndarray::tables::Index<3> const PARAMETERS;
    static ndarray::tables::Index<4> const COEFFICIENTS;
    static ndarray::tables::Index<4> const FISHER_INFO;

    typedef ndarray::tables::Table<ImportanceSample> Table;
    typedef ndarray::tables::Layout<ImportanceSample> Layout;
};

}}}} // namespace lsst::meas::multifit::sampling

#endif // !LSST_MEAS_MULTIFIT_SAMPLING_ImportanceSample
