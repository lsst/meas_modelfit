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
 
#ifndef LSST_MEAS_MULTIFIT_CACHE_H
#define LSST_MEAS_MULTIFIT_CACHE_H

#include <Eigen/Core>

#include "lsst/afw/math/Function.h"
#include "lsst/afw/geom.h"

namespace boost { 
namespace serialization {
    class access;
}}

namespace lsst{
namespace meas{
namespace multifit {

#ifndef SWIG
class InterpolationFunction {
public:
    typedef boost::shared_ptr<InterpolationFunction> Ptr;
    typedef boost::shared_ptr<InterpolationFunction const> ConstPtr;

    InterpolationFunction(
        double const min, double const max, double const step,
        Eigen::VectorXd const & values 
    ) : _min(min), _max(max), _step(step), _values(values) {}

    virtual double operator() (double x) const;
    virtual double d(double x) const;

protected:
    double _min, _max, _step;
    Eigen::VectorXd _values;
};

class InterpolationFunctionFactory : private boost::noncopyable {
public:
    
    typedef boost::shared_ptr<InterpolationFunctionFactory> Ptr;
    typedef boost::shared_ptr<InterpolationFunctionFactory const> ConstPtr;

    struct Registration {
        explicit Registration(Ptr const & p);
    };

    static InterpolationFunctionFactory::ConstPtr get(std::string const & name);

    std::string const & getName() const { return _name; }

    int getExtraParameterSize() const { return _extraParameterSize; }
    
    InterpolationFunction::ConstPtr operator()(
        double const min, double const max, double const step,
        Eigen::VectorXd const & values
    ) const;

    virtual void fillExtraParameters(
        double const fastMax, 
        double const slow,
        Eigen::MatrixXd::RowXpr row, 
        lsst::afw::math::Function2<double> const & fillFunction
    ) const {}

    virtual ~InterpolationFunctionFactory() {}

protected:
    explicit InterpolationFunctionFactory(std::string const & name, int extraParameterSize);

    virtual InterpolationFunction::ConstPtr execute(
        double const min, double const mac, double const step,
        Eigen::VectorXd const & values
    ) const = 0;

    typedef std::map<std::string, InterpolationFunctionFactory::Ptr> Registry; 

    std::string _name;
    int _extraParameterSize;

private:
    static Registry & getRegistry();
};

#endif

class Cache  {
public:   
    typedef boost::shared_ptr<Cache> Ptr;
    typedef boost::shared_ptr<Cache const> ConstPtr;
    typedef lsst::afw::math::Function2<double> FillFunction;

    typedef InterpolationFunction Interpolator;
    typedef InterpolationFunctionFactory InterpolatorFactory;

    static ConstPtr get(std::string const & name);
    static ConstPtr make(
        double const slowMin, double const slowMax,
        double const fastMin, double const fastMax,
        int const nSlow, int const nFast,
        FillFunction const * fillFunction,
        InterpolatorFactory::ConstPtr const & interpolatorFactory,
        std::string const & name="", 
        bool const & doOverwrite=true
    );
    static ConstPtr load(
        std::string const & filepath, 
        std::string const & name="", 
        bool const & doOverwrite=true
    ); 
    void save(std::string const & filepath) const;

#ifndef SWIG
    Cache::Interpolator::ConstPtr getInterpolator(double const slow) const;
    Cache::Interpolator::ConstPtr getDerivativeInterpolator(double const slow) const;
#endif

    double const getSlowMin() const {return _slowMin;}
    double const getSlowMax() const {return _slowMax;}
    double const getFastMin() const {return _fastMin;}
    double const getFastMax() const {return _fastMax;}
    double const getSlowStep() const {return _slowStep;}
    double const getFastStep() const {return _fastStep;}
    int const getNSlow() const {return int((_slowMax - _slowMin)/_slowStep);}
    int const getNFast() const {return int((_fastMax - _fastMin)/_fastStep);}

    Eigen::MatrixXd const & getDataPoints() const {return _dataPoints;}

private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & ar, unsigned int const);
    
    Cache(
        double const slowMin, double const slowMax,
        double const fastMin, double const fastMax,
        int const nSlow, int const nFast,
        FillFunction const * fillFunction,
        InterpolatorFactory::ConstPtr const & interpolatorFactory
    );

    Cache(){};

    Eigen::MatrixXd _dataPoints;

    InterpolatorFactory::ConstPtr _interpolatorFactory;

    double _slowMin, _slowMax, _fastMin, _fastMax;
    double _slowStep, _fastStep;

    static std::map<std::string, Ptr> _registry;
};

}}} //namespace lsst::meas::multifit
#endif
