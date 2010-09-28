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
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & y
    ) : Base(x, y) {}

    RationalExtrapolationInterpolator (
        Eigen::VectorXd const & x,
        std::vector<double> const & y
    ) : Base(x, y) {}

    virtual Base::Base::Ptr clone() const {
        return boost::make_shared<RationalExtrapolationInterpolator>(_x, _params);
    }

    virtual double operator()(double x) const {
        if (x < _x[_x.size() - 1]) return Base::operator()(x);
        return _params[_x.size()] + _params[_x.size() + 1] / x + _params[_x.size() + 2] / (x*x);
    }

    virtual double d(double x) const {
        if (x < _x[_x.size() - 1]) return Base::d(x);
        return -_params[_x.size() + 1] / (x*x) - 2.0 * _params[_x.size() + 2] / (x*x*x);
    }

};

class RationalExtrapolationInterpolatorFactory : public multifit::InterpolationFunctionFactory {
public:

    static Registration registration;

    virtual void fillExtraParameters(
        Eigen::VectorXd const & x, 
        double y,
        Eigen::MatrixXd::RowXpr row, 
        lsst::afw::math::Function2<double> const & fillFunction
    ) const {
        if (!_complete) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Cannot fill cache with incomplete FunctorFactory."
            );
        }
        double x0 = x[x.size() - 1];
        double z0 = row[x.size() - 1];
        double z1 = fillFunction(_x1, y);
        double z2 = fillFunction(_x2, y);
        fitRational(
            Eigen::Vector3d(x0, _x1, _x2),
            Eigen::Vector3d(z0, z1, z2),
            row[row.size() - 3], row[row.size() - 2], row[row.size() - 1]
        );
    }

    virtual void fillExtraParameters(
        double x,
        Eigen::VectorXd const & y, 
        Eigen::MatrixXd::ColXpr col, 
        lsst::afw::math::Function2<double> const & fillFunction
    ) const {
        if (!_complete) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Cannot fill cache with incomplete FunctorFactory."
            );
        }
        double x0 = y[y.size() - 1];
        double z0 = col[y.size() - 1];
        double z1 = fillFunction(x, _x1);
        double z2 = fillFunction(x, _x2);
        fitRational(
            Eigen::Vector3d(x0, _x1, _x2),
            Eigen::Vector3d(z0, z1, z2),
            col[col.size() - 3], col[col.size() - 2], col[col.size() - 1]
        );
    }

    RationalExtrapolationInterpolatorFactory() :
        multifit::InterpolationFunctionFactory("RationalExtrapolation", 3),
        _complete(false), _x1(0.0), _x2(0.0)
    {}

    RationalExtrapolationInterpolatorFactory(double x1, double x2) :
        multifit::InterpolationFunctionFactory("RationalExtrapolation", 3),
        _complete(true), _x1(x1), _x2(x2)
    {}

protected:

    virtual multifit::InterpolationFunction::ConstPtr execute(
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & y
    ) const {
        return boost::make_shared<RationalExtrapolationInterpolator>(x, y);
    }

private:

    static void fitRational(
        Eigen::Vector3d const & x, Eigen::Vector3d const & y,
        double & a0, double & a1, double & a2
    ) {
        Eigen::Matrix3d matrix;
        matrix.col(0).fill(1.0);
        matrix.col(1) = x.cwise().inverse();
        matrix.col(2) = x.cwise().inverse().cwise().square();
        Eigen::Vector3d model;
        Eigen::LU<Eigen::Matrix3d> lu(matrix);
        lu.solve(y, &model);
        a0 = model[0];
        a1 = model[1];
        a2 = model[2];
        if (lsst::utils::isnan(a0) || lsst::utils::isnan(a1) || lsst::utils::isnan(a2))
            a0 = a1 = a2 = 0.0;        
    }

    bool _complete;
    double _x1, _x2;
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
    Cache::FunctorFactory::Ptr rowFunctorFactory
        = boost::make_shared<RationalExtrapolationInterpolatorFactory>(
            policy.getDouble("extrapolationK1"),
            policy.getDouble("extrapolationK2")
        );
    RobustSersicCacheFillFunction fillFunction(
        policy.getDouble("sersicInnerRadius"), 
        policy.getDouble("sersicOuterRadius"), 
        policy.getDouble("epsabs"), 
        policy.getDouble("epsrel"),
        policy.getBool("noInterpolation")
    );

    return Cache::make(
        bounds, resolution, &fillFunction, 
        rowFunctorFactory, Cache::FunctorFactory::get(""),
        "RobustSersic", false
    );
}
