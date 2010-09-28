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
class InterpolationFunction : public lsst::afw::math::Function1<double> {
public:
    typedef lsst::afw::math::Function1<double> Base;
    typedef boost::shared_ptr<InterpolationFunction> Ptr;
    typedef boost::shared_ptr<InterpolationFunction const> ConstPtr;

    InterpolationFunction(
        Eigen::VectorXd const & x, 
        std::vector<double> const & y
    );
    InterpolationFunction(
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & y 
    );

    virtual Base::Ptr clone() const;
    virtual double operator() (double x) const;
    virtual double d(double x) const;

protected:

    Eigen::VectorXd _x;
    double _step;
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
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & y
    ) const;

    virtual void fillExtraParameters(
        Eigen::VectorXd const & x, 
        double y,
        Eigen::MatrixXd::RowXpr row, 
        lsst::afw::math::Function2<double> const & fillFunction
    ) const {}

    virtual void fillExtraParameters(
        double x,
        Eigen::VectorXd const & y, 
        Eigen::MatrixXd::ColXpr col, 
        lsst::afw::math::Function2<double> const & fillFunction
    ) const {}

    virtual ~InterpolationFunctionFactory() {}

protected:
    explicit InterpolationFunctionFactory(std::string const & name, int extraParameterSize);

    virtual InterpolationFunction::ConstPtr execute(
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & y
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

    typedef InterpolationFunction Functor;
    typedef InterpolationFunctionFactory FunctorFactory;

    static ConstPtr get(std::string const & name);
    static ConstPtr make( 
        lsst::afw::geom::BoxD const & parameterBounds,
        lsst::afw::geom::Extent2D const & resolution,
        FillFunction const * fillFunction,
        FunctorFactory::ConstPtr const & rowFunctorFactory,
        FunctorFactory::ConstPtr const & colFunctorFactory,
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
    Cache::Functor::ConstPtr getRowFunctor(double const & y) const;
    Cache::Functor::ConstPtr getRowDerivativeFunctor(double const & y) const;
    Cache::Functor::ConstPtr getColFunctor(double const & x) const;
    Cache::Functor::ConstPtr getColDerivativeFunctor(double const & x) const;

    Eigen::VectorXd const & getColHeaders() const {return _x;}
    Eigen::VectorXd const & getRowHeaders() const {return _y;}
    
    lsst::afw::geom::BoxD const & getParameterBounds() const {
        return _parameterBounds;
    }
#endif
    Eigen::MatrixXd const & getDataPoints() const {return _dataPoints;}
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & ar, unsigned int const);
    
    Cache(
        lsst::afw::geom::BoxD const & parameterBounds,
        lsst::afw::geom::Extent2D const & resolution,
        FillFunction const * fillFunction,
        FunctorFactory::ConstPtr const & rowFunctorFactory,
        FunctorFactory::ConstPtr const & colFunctorFactory
    );

    Cache(){};

    Eigen::MatrixXd _dataPoints;
    Eigen::VectorXd _x, _y;

    lsst::afw::geom::BoxD _parameterBounds;
    double _xStep, _yStep;
    
    FunctorFactory::ConstPtr _rowFunctorFactory;
    FunctorFactory::ConstPtr _colFunctorFactory;

    static std::map<std::string, Ptr> _registry;
};

}}} //namespace lsst::meas::multifit
#endif
