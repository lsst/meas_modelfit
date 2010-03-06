#ifndef LSST_MEAS_MULTIFIT_CACHE_H
#define LSST_MEAS_MULTIFIT_CACHE_H

#include <Eigen/Core>

#include "lsst/afw/math/Function.h"
#include "lsst/afw/geom.h"

namespace lsst{
namespace meas{
namespace multifit {


class Cache :public lsst::daf::base::Persistable {
public:   
    typedef boost::shared_ptr<Cache> Ptr;
    typedef boost::shared_ptr<Cache const> ConstPtr;
    typedef lsst::afw::math::Function2<double> FillFunction;

    class Functor : public lsst::afw::math::Function1<double> {
    public:
        typedef lsst::afw::math::Function1<double> Base;
        typedef boost::shared_ptr<Functor> Ptr;
        typedef boost::shared_ptr<Functor const> ConstPtr;

        Functor(
            boost::shared_ptr<Eigen::VectorXd> const & x, 
            std::vector<double> const & y
        );
        Functor(
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

    Cache(
        lsst::afw::geom::BoxD const & parameterBounds,
        lsst::afw::geom::Extent2D const & resolution,
        FillFunction::Ptr const & fillFunction
    );

    Functor::ConstPtr getRowFunctor(double const & y) const;
    Functor::ConstPtr getRowDerivativeFunctor(double const & y) const;
    Functor::ConstPtr getColFunctor(double const & x) const;
    Functor::ConstPtr getColDerivativeFunctor(double const & x) const;

    Eigen::VectorXd const & getColHeaders() const {return *_x;}
    Eigen::VectorXd const & getRowHeaders() const {return *_y;}
    
    lsst::afw::geom::BoxD const & getParameterBounds() const {
        return _parameterBounds;
    }
private:
    Eigen::MatrixXd _dataPoints;
    boost::shared_ptr<Eigen::VectorXd> _x, _y;

    lsst::afw::geom::BoxD _parameterBounds;
    double _xStep, _yStep;
    static double const _getStep(
        int const &i, 
        boost::shared_ptr<Eigen::VectorXd> r
    ) {
        return (*r)[i+1] - (*r)[i];
    }
};

}}} //namespace lsst::meas::multifit
#endif
