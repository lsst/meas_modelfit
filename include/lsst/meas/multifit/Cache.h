#ifndef LSST_MEAS_MULTIFIT_CACHE_H
#define LSST_MEAS_MULTIFIT_CACHE_H

#include <Eigen/Core>

#include "lsst/afw/math/Function.h"
#include "lsst/afw/geom.h"

namespace lsst{
namespace meas{
namespace multifit {

class InterpolationFunction : public lsst::afw::math::Function1<double> {
public:
    typedef lsst::afw::math::Function1<double> Base;
    typedef boost::shared_ptr<InterpolationFunction> Ptr;
    typedef boost::shared_ptr<InterpolationFunction const> ConstPtr;

    InterpolationFunction(
        boost::shared_ptr<Eigen::VectorXd> const & x, 
        std::vector<double> const & y
    );
    InterpolationFunction(
        boost::shared_ptr<Eigen::VectorXd> const & x,
        Eigen::VectorXd const & y 
    );

    virtual Base::Ptr clone() const;
    virtual double operator() (double x) const;
    virtual double dParams(double const & x) const;

private:
    void initialize();

    boost::shared_ptr<Eigen::VectorXd> _x;
    double _step;
}; 

class Cache {
public:    
    typedef lsst::afw::math::Function2<double> FillFunction;

    Cache(
        lsst::afw::geom::Extent2I const & dimensions,
        lsst::afw::geom::BoxD const & parameterBounds,
        FillFunction::Ptr const & fillFunction
    );

    InterpolationFunction::ConstPtr getRowFunctor(double const & y) const;
    InterpolationFunction::ConstPtr getRowDerivativeFunctor(double const & y) const;
    InterpolationFunction::ConstPtr getColFunctor(double const & x) const;
    InterpolationFunction::ConstPtr getColDerivativeFunctor(double const & x) const;

    Eigen::VectorXd const & getColHeaders() const {return *_x;}
    Eigen::VectorXd const & getRowHeaders() const {return *_y;}
    
    lsst::afw::geom::BoxD const & getParameterBounds() const {
        return _parameterBounds;
    }
    double const & getColStep() const {return _xStep;}
    double const & getRowStep() const {return _yStep;}

private:    
    Eigen::MatrixXd _dataPoints;
    boost::shared_ptr<Eigen::VectorXd> _x, _y;

    lsst::afw::geom::BoxD _parameterBounds;
    double _xStep, _yStep;
};

}}} //namespace lsst::meas::multifit
#endif
