#include "lsst/meas/multifit/Cache.h"

#include <cmath>
#include <boost/make_shared.hpp>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::Cache::Cache(
    Eigen::MatrixXd const & dataPoints,
    lsst::afw::geom::BoxD parameterBounds
) : _dataPoints(dataPoints),
    _x(new Eigen::VectorXd(dataPoints.cols())), 
    _y(new Eigen::VectorXd(dataPoints.rows())),
    _parameterBounds(parameterBounds),
    _xStep(parameterBounds.getWidth() / dataPoints.cols()),
    _yStep(parameterBounds.getHeight() / dataPoints.rows())
{
    if(dataPoints.rows() < 2 || dataPoints.cols() < 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "At least 2 data points needed in each parameter axis"
        );
    }

    double x = parameterBounds.getMinX();
    double y = parameterBounds.getMaxX();
    for(int i = 0; i < _x->size() && i < _y->size(); ++i) {
        (*_x)[i] = x;
        (*_y)[i] = y;
        x += _xStep;
        y += _yStep;
    }
    //if more cols than rows
    //fill in remaining column headers
    for(int i = _y->size(); i < _x->size(); ++i) {
        (*_x)[i] = x;
        x += _xStep;
    }
    //if more rows than cols
    //fill in remining row headers
    for(int i = _x->size(); i < _y->size(); ++i) {
        (*_y)[i] = y;
        y += _yStep;
    }
};

multifit::Cache::Cache(
    Eigen::MatrixXd const & dataPoints,
    Eigen::VectorXd const & x,
    Eigen::VectorXd const & y
) : _dataPoints(dataPoints), 
    _x(new Eigen::VectorXd(x)), 
    _y(new Eigen::VectorXd(y)) 
{
    if(dataPoints.rows() < 2 || dataPoints.cols() < 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "At least 2 data points needed in each parameter axis"
        );
    }
    if(_dataPoints.cols() != _x->size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Length of x must match number of columns in dataPoints"
        );
    }  
    if(_dataPoints.rows() != _y->size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Length of y must match number of rows in dataPoints"
        );
    }  
   
    lsst::afw::geom::PointD min = lsst::afw::geom::makePointD(
        (*_x)[0], (*_y)[0]
    );
    lsst::afw::geom::PointD max = lsst::afw::geom::makePointD(
        (*_x)[_x->size() -1], (*_y)[_y->size() - 1]
    );
    _parameterBounds = lsst::afw::geom::BoxD(min, max);
    _xStep = _parameterBounds.getWidth() / _dataPoints.cols();
    _yStep = _parameterBounds.getHeight() / _dataPoints.rows();
} 

multifit::Cache::Cache(
    lsst::afw::geom::Extent2I const & dimensions,
    lsst::afw::geom::BoxD const & parameterBounds,
    bool computeHeaders
) : _dataPoints(dimensions.getY(), dimensions.getX()),
    _x(new Eigen::VectorXd(dimensions.getX())),
    _y(new Eigen::VectorXd(dimensions.getY())),
    _parameterBounds(parameterBounds),
    _xStep(parameterBounds.getWidth()/dimensions.getX()),
    _yStep(parameterBounds.getHeight()/dimensions.getY())
{
    if(dimensions.getX() < 2 || dimensions.getY() < 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "At least 2 data points needed in each parameter axis"
        );
    }

    if(computeHeaders) { 
        double x = _parameterBounds.getMinX();
        double y = _parameterBounds.getMinY();
        for(int i = 0; i < _x->size() && i < _y->size(); ++i) {
            (*_x)[i] = x;
            (*_y)[i] = y;
            x += _xStep;
            y += _yStep;
        }
        //if more cols than rows
        //fill in remaining column headers
        for(int i = _y->size(); i < _x->size(); ++i) {
            (*_x)[i] = x;
            x += _xStep;
        }
        //if more rows than cols
        //fill in remining row headers
        for(int i = _x->size(); i < _y->size(); ++i) {
            (*_y)[i] = y;
            y += _yStep;
        } 
    }
}


multifit::InterpolationFunction::ConstPtr multifit::Cache::getRowFunctor(
    double const & y
) const {
    int i = static_cast<int>(y / _yStep);
    double s = (y - (*_y)[i]) / _yStep;
   
    Eigen::VectorXd r(_dataPoints.row(i)*(1-s) + _dataPoints.row(i+1)*s);
    return boost::make_shared<InterpolationFunction const>(_y, r);
}

multifit::InterpolationFunction::ConstPtr multifit::Cache::getRowDerivativeFunctor(
    double const & y
) const {
    int i = static_cast<int>(y / _yStep);

    Eigen::VectorXd r((_dataPoints.row(i) - _dataPoints.row(i))/_yStep);
    return boost::make_shared<InterpolationFunction const>(_y, r);
}

multifit::InterpolationFunction::ConstPtr multifit::Cache::getColFunctor(
    double const & x
) const {
    int i = static_cast<int>(x / _xStep);
    double s = (x - (*_x)[i]) / _xStep;

    Eigen::VectorXd r(_dataPoints.col(i)*(1-s) + _dataPoints.col(i+1)*s);
    return boost::make_shared<InterpolationFunction const>(_x, r);
}

multifit::InterpolationFunction::ConstPtr multifit::Cache::getColDerivativeFunctor(
    double const & x
) const {
    int i = static_cast<int>(x / _xStep);

    Eigen::VectorXd r((_dataPoints.col(i) - _dataPoints.col(i))/_xStep);
    return boost::make_shared<InterpolationFunction const>(_x, r);
}

double multifit::InterpolationFunction::operator() (double x) const {
    int i = static_cast<int>(x / _step);
    double w1 = (x - (*_x)[i])/_step;
    double w0 = 1 - w1;
    double p0 = _params[i];
    double p1 = _params[i+1];

    return w0*p0+w1*p1;
}

double multifit::InterpolationFunction::dParams(double const & x) const {
    int i = static_cast<int>(x / _step);
    double p0 = _params[i];
    double p1 = _params[i+1];

    return (p1 - p0)/ _step;
}

multifit::InterpolationFunction::InterpolationFunction(
    boost::shared_ptr<Eigen::VectorXd> const & x, 
    std::vector<double> const & y
) : Base(y), _x(x) {
    initialize();  
}

multifit::InterpolationFunction::InterpolationFunction(
    boost::shared_ptr<Eigen::VectorXd> const & x,
    Eigen::VectorXd const & y 
) : Base(y.size()), _x(x) {
    initialize();
    std::copy(y.data(), y.data() + y.size(), _params.begin());
}

multifit::InterpolationFunction::Base::Ptr multifit::InterpolationFunction::clone() const {
    return boost::make_shared<InterpolationFunction>(_x, _params);
}
void multifit::InterpolationFunction::initialize() {
    if(_x->size() < 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Not enough points. At least 2 data points needed"
        );
    }
    if(getNParameters() != static_cast<unsigned int>(_x->size())) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Length of vector x must match length of vector y"
        );
    }

    _step = (*_x)[1] - (*_x)[0];
}

