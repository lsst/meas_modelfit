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
    virtual double dParams(double const & x) const;

private:
    void initialize();

    Eigen::VectorXd _x;
    double _step;
}; 
#endif

class Cache  {
public:   
    typedef boost::shared_ptr<Cache> Ptr;
    typedef boost::shared_ptr<Cache const> ConstPtr;
    typedef lsst::afw::math::Function2<double> FillFunction;

    typedef InterpolationFunction Functor;

    static ConstPtr get(std::string const & name);
    static ConstPtr make( 
        lsst::afw::geom::BoxD const & parameterBounds,
        lsst::afw::geom::Extent2D const & resolution,
        FillFunction const * fillFunction,
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

    Eigen::MatrixXd _dataPoints;
    Eigen::VectorXd _x, _y;
    
    Cache(
        lsst::afw::geom::BoxD const & parameterBounds,
        lsst::afw::geom::Extent2D const & resolution,
        FillFunction const * fillFunction
    );

    Cache(){};

    lsst::afw::geom::BoxD _parameterBounds;
    double _xStep, _yStep;
    


    static std::map<std::string, Ptr> _registry;
};

}}} //namespace lsst::meas::multifit
#endif
