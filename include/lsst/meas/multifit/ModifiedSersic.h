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
 
#ifndef LSST_MEAS_MULTIFIT_MODIFIED_SERSIC_H
#define LSST_MEAS_MULTIFIT_MODIFIED_SERSIC_H

#include <boost/make_shared.hpp>
#include "lsst/afw/math/Function.h"
#include "lsst/pex/policy/DefaultPolicyFile.h"

namespace lsst {
namespace meas {
namespace multifit {

class ModifiedSersicFunction {
public:

    ModifiedSersicFunction(double n, double inner, double outer);

    double operator()(double r) const;

    double getSersicIndex() const { return  _n; }
    void setSersicIndex(double n);
    double getOuter() const { return _outer; }
    double getInner() const { return _inner; }

    double integrate(double radius, int moment) const;

private:

    double integrateInner(double radius, int moment) const;
    double integrateOuter(double radius, int moment) const;

    double _n;
    double _inner, _outer;
    double _a[4];
};

class ModifiedSersicHankelTransform {
public:

    ModifiedSersicHankelTransform(
        ModifiedSersicFunction const & function, 
        double epsrel=-1.0,  ///< relative tolerance; defaults to sqrt(epsilon) if < 0.0
        double epsabs=-1.0,  ///< absolute tolerance; defaults to epsilon if < 0.0,
        bool doInterpolation=false
    );

    double operator()(double k) const;

    ModifiedSersicFunction const & getFunction() const { return  _function; }
    
    void setFunction(ModifiedSersicFunction const & function) { _function = function; }
    void setSersicIndex(double n) {_function.setSersicIndex(n); }
private:

    /// Robust sum of a vector by sorting first.  Will not preserve ordering of the elements!
    static double sum(std::vector<double> & v);

    ModifiedSersicFunction _function;
    double _epsrel;
    double _epsabs;
    bool _doInterpolation;
    mutable std::vector<double> _endpoints;
    mutable std::vector<double> _integrals;
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_MODIFIED_SERSIC_H
