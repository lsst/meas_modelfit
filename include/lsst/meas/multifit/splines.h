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

#ifndef LSST_MEAS_MULTIFIT_splines_h_INCLUDED
#define LSST_MEAS_MULTIFIT_splines_h_INCLUDED

#include "lsst/base.h"
#include "ndarray.h"

namespace lsst { namespace meas { namespace multifit {

class BasisSpline {
public:

    explicit BasisSpline(ndarray::Array<double const,1,1> const & knots, int degree=3);

    int getBasisSize() const;

    void evaluate(double x, ndarray::Array<double,1,1> const & b) const;

    ndarray::Array<double,1,1> evaluate(double x) const;

    void evaluate(
        ndarray::Array<double const,1,1> const & x,
        ndarray::Array<double,2,1> const & b
    ) const;

    ndarray::Array<double,2,2> evaluate(
        ndarray::Array<double const,1,1> const & x
    ) const;

    void evaluateDerivatives(double x, ndarray::Array<double,2,1> const & db) const;

    ndarray::Array<double,2,2> evaluateDerivatives(double x, int nDeriv) const;

    ~BasisSpline(); // needs to be in .cc so compiler can see ~Impl()

private:
    class Impl;
    PTR(Impl) _impl;
};


}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_splines_h_INCLUDED
