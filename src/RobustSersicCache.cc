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
 
#include "lsst/meas/multifit/RobustSersicCache.h"

namespace multifit = lsst::meas::multifit;

multifit::RobustSersicCacheFillFunction::RobustSersicCacheFillFunction(
    double inner, double outer, double epsabs, double epsrel, bool noInterpolation
) : 
    Cache::FillFunction(0),
    _hankelTransform(ModifiedSersicFunction(1.0, inner, outer), epsrel, epsabs, noInterpolation),
    _epsrel(epsrel), _epsabs(epsabs), _noInterpolation(noInterpolation)
{}

double multifit::RobustSersicCacheFillFunction::operator() (double x, double y) const {
    if(y != _hankelTransform.getFunction().getSersicIndex()) {
        _hankelTransform.setFunction(
            ModifiedSersicFunction(
                y, 
                _hankelTransform.getFunction().getInner(), 
                _hankelTransform.getFunction().getOuter()
            )
        );
    } 
    return _hankelTransform(x);
}

multifit::Cache::ConstPtr multifit::makeRobustSersicCache(lsst::pex::policy::Policy policy) {
    lsst::pex::policy::DefaultPolicyFile defSource(
        "meas_multifit",
        "RobustSersicCacheDict.paf",
        "policy"
    );
    lsst::pex::policy::Policy defPol(defSource);
    policy.mergeDefaults(defPol);

    lsst::afw::geom::Extent2D resolution = lsst::afw::geom::makeExtentD(
        policy.getDouble("kResolution"), 
        policy.getDouble("sersicIndexResolution")
    );
    lsst::afw::geom::BoxD bounds(
        lsst::afw::geom::makePointD(
            policy.getDouble("kMin"), 
            policy.getDouble("sersicIndexMin")
        ),
        lsst::afw::geom::makePointD(
            policy.getDouble("kMax"),
            policy.getDouble("sersicIndexMax")
        )
    );
    RobustSersicCacheFillFunction fillFunction(
        policy.getDouble("inner"), 
        policy.getDouble("outer"), 
        policy.getDouble("epsabs"), 
        policy.getDouble("epsrel"),
        policy.getBool("noInterpolation")
    );

    return Cache::make(bounds, resolution, &fillFunction, "", "", "RobustSersic", false);
}
