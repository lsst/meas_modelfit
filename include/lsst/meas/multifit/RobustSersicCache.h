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
 
#ifndef LSST_MEAS_MULTIFIT_ROBUST_SERSIC_CACHE_FILL_FUNCTION_H
#define LSST_MEAS_MULTIFIT_ROBUST_SERSIC_CACHE_FILL_FUNCTION_H

#include <boost/make_shared.hpp>
#include "lsst/afw/geom.h"
#include "lsst/meas/multifit/Cache.h"
#include "lsst/meas/multifit/ModifiedSersic.h"
#include "lsst/pex/policy/DefaultPolicyFile.h"

namespace lsst {
namespace meas {
namespace multifit {

#ifndef SWIG
class RobustSersicCacheFillFunction : public Cache::FillFunction {
public: 

    RobustSersicCacheFillFunction(
        double inner, double outer, 
        double epsrel=-1.0, double epsabs=-1.0,
        bool noInterpolation=false
    );

    virtual double operator()(double x, double y) const;

    virtual Cache::FillFunction::Ptr clone() const {
        return boost::make_shared<RobustSersicCacheFillFunction>(
            _hankelTransform.getFunction().getInner(), _hankelTransform.getFunction().getOuter(), 
            _epsrel, _epsabs, _noInterpolation
        );
    }

private:
    mutable ModifiedSersicHankelTransform _hankelTransform;
    double _epsrel;
    double _epsabs;
    bool _noInterpolation;
};
#endif

Cache::ConstPtr makeRobustSersicCache(lsst::pex::policy::Policy policy);

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_ROBUST_SERSIC_CACHE_FILL_FUNCTION_H
