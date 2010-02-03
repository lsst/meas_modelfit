#include "lsst/meas/multifit/Cache.h"

#include <cmath>
#include <boost/make_shared.hpp>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::Cache::Cache(
    lsst::afw::geom::Extent2I const& dimensions,
    lsst::afw::geom::BoxD const & parameterBounds,
    FillFunction::Ptr const & fillFunction
) : _dataPoints(dimensions.getY(), dimensions.getX()),
    _x(new Eigen::VectorXd(dimensions.getX())), 
    _y(new Eigen::VectorXd(dimensions.getY())),
    _parameterBounds(parameterBounds),
    _xStep(parameterBounds.getWidth() / dimensions.getX()),
    _yStep(parameterBounds.getHeight() / dimensions.getY())
{
    if(dimensions.getX() < 2 || dimensions.getY() < 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "At least 2 data points needed in each parameter axis"
        );
    }

    //loop over columns to precompute headers
    double x = _parameterBounds.getMinX();
    double * xHeader = _x->data();
    for(int i = 0; i < dimensions.getX(); ++i, ++xHeader) {
        *xHeader = x;
        x += _xStep;
    } 

    double * yHeader = _y->data();
    double y = _parameterBounds.getMinY();
    
    //outer loop over rows
    for(int i = 0; i < dimensions.getY(); ++i, ++yHeader) {
        *yHeader = y;
        
        x = _parameterBounds.getMinX();
        //inner loop over columns
        for (int j = 0; j < dimensions.getX(); ++j) {
            //set grid point
            _dataPoints(i,j) = (*fillFunction)(x, y);

            x += _xStep;
        }

        y += _yStep;
    }
};

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

