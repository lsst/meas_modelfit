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
#include "lsst/pex/exceptions.h"
#include "lsst/utils/ieee.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Array>

namespace multifit = lsst::meas::multifit;

namespace {

class RationalExtrapolationInterpolator : public multifit::InterpolationFunction {
public:

    typedef multifit::InterpolationFunction Base;
    
    RationalExtrapolationInterpolator (
        double const min, double const max, double const step,
        Eigen::VectorXd const & values
    ) : Base(min, max, step, values) {}

    virtual double operator()(double p) const {
        if (p < _max) return Base::operator()(p);
        int n = _values.size();
        return _values[n-3] + _values[n-2] / p + _values[n-1] / (p*p);
    }

    virtual double d(double p) const {
        if (p < _max) return Base::d(p);
        int n = _values.size();
        return -_values[n-2] / (p*p) - 2.0 * _values[n-1] / (p*p*p);
    }

};

class RationalExtrapolationInterpolatorFactory : public multifit::InterpolationFunctionFactory {
public:

    static Registration registration;

    virtual void fillExtraParameters(
        double fastMax, 
        double slow,
        Eigen::MatrixXd::RowXpr row, 
        lsst::afw::math::Function2<double> const & fillFunction
    ) const {
        if (!_complete) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Cannot fill cache with incomplete InterpolatorFactory."
            );
        }
        double v0 = row[row.size() - getExtraParameterSize() - 1];
        double v1 = fillFunction(slow, _fast1);
        double v2 = fillFunction(slow, _fast2);
        fitRational(
            Eigen::Vector3d(fastMax, _fast1, _fast2),
            Eigen::Vector3d(v0, v1, v2),
            row[row.size() - 3], row[row.size() - 2], row[row.size() - 1]
        );
    }

    RationalExtrapolationInterpolatorFactory() :
        multifit::InterpolationFunctionFactory("RationalExtrapolation", 3),
        _complete(false), _fast1(0.0), _fast2(0.0)
    {}

    RationalExtrapolationInterpolatorFactory(double fast1, double fast2) :
        multifit::InterpolationFunctionFactory("RationalExtrapolation", 3),
        _complete(true), _fast1(fast1), _fast2(fast2)
    {}

protected:

    virtual multifit::InterpolationFunction::ConstPtr execute(
        double const min, double const max, double const step,
        Eigen::VectorXd const & values
    ) const {
        return boost::make_shared<RationalExtrapolationInterpolator>(
            min, max, step, values
        );
    }

private:

    static void fitRational(
        Eigen::Vector3d const & fast, Eigen::Vector3d const & values,
        double & a0, double & a1, double & a2
    ) {
        Eigen::Matrix3d matrix;
        matrix.col(0).fill(1.0);
        matrix.col(1) = fast.cwise().inverse();
        matrix.col(2) = fast.cwise().inverse().cwise().square();
        Eigen::Vector3d model;
        Eigen::LU<Eigen::Matrix3d> lu(matrix);
        lu.solve(values, &model);
        a0 = model[0];
        a1 = model[1];
        a2 = model[2];
        if (lsst::utils::isnan(a0) || lsst::utils::isnan(a1) || lsst::utils::isnan(a2))
            a0 = a1 = a2 = 0.0;        
    }

    bool _complete;
    double _fast1, _fast2;
};

RationalExtrapolationInterpolatorFactory::Registration
RationalExtrapolationInterpolatorFactory::registration(
    boost::make_shared<RationalExtrapolationInterpolatorFactory>()
);

} // anonymous

multifit::RobustSersicCacheFillFunction::RobustSersicCacheFillFunction(
    double inner, double outer, double epsabs, double epsrel, bool noInterpolation
) : 
    Cache::FillFunction(0),
    _hankelTransform(ModifiedSersicFunction(1.0, inner, outer), epsrel, epsabs, noInterpolation),
    _epsrel(epsrel), _epsabs(epsabs), _noInterpolation(noInterpolation)
{}

double multifit::RobustSersicCacheFillFunction::operator() (double sersic, double k) const {
    if(sersic != _hankelTransform.getFunction().getSersicIndex()) {
        _hankelTransform.setFunction(
            ModifiedSersicFunction(
                sersic, 
                _hankelTransform.getFunction().getInner(), 
                _hankelTransform.getFunction().getOuter()
            )
        );
    } 
    return _hankelTransform(k);
}

multifit::Cache::ConstPtr multifit::makeRobustSersicCache(lsst::pex::policy::Policy policy) {
    lsst::pex::policy::DefaultPolicyFile defSource(
        "meas_multifit",
        "RobustSersicCacheDict.paf",
        "policy"
    );
    lsst::pex::policy::Policy defPol(defSource);
    policy.mergeDefaults(defPol);

    Cache::InterpolatorFactory::Ptr interpolatorFactory(
        new RationalExtrapolationInterpolatorFactory(
            policy.getDouble("extrapolationK1"),
            policy.getDouble("extrapolationK2")
        )
    );
    RobustSersicCacheFillFunction fillFunction(
        policy.getDouble("sersicInnerRadius"), 
        policy.getDouble("sersicOuterRadius"), 
        policy.getDouble("epsabs"), 
        policy.getDouble("epsrel"),
        policy.getBool("noInterpolation")
    );

    return Cache::make(
        policy.getDouble("sersicMin"), policy.getDouble("sersicMax"),
        policy.getDouble("kMin"), policy.getDouble("kMax"),
        policy.getInt("nBinSersic"), policy.getInt("nBinK"),
        &fillFunction, 
        interpolatorFactory, 
        "RobustSersic", 
        false
    );
}
