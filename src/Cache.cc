#include "lsst/meas/multifit/Cache.h"

#include <cmath>
#include <boost/make_shared.hpp>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::Cache::Cache(
    lsst::afw::geom::BoxD const & parameterBounds,
    lsst::afw::geom::Extent2D const & resolution,
    FillFunction::Ptr const & fillFunction
) : _parameterBounds(parameterBounds),
    _xStep(resolution.getX()),
    _yStep(resolution.getY())
{
    if(_xStep > parameterBounds.getWidth() || 
       _yStep > parameterBounds.getHeight()
    ) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Requested resolution exceed size of parameter Bounds"
        );
    }
    
    int nRow = static_cast<int>(_parameterBounds.getHeight()/_yStep);
    int nCol = static_cast<int>(_parameterBounds.getWidth()/_xStep);
    
    //allocate necesary space (note the +1 in order to be inclusive of upper
    //parameter bounds
    _dataPoints = Eigen::MatrixXd(nRow + 1, nCol + 1);
    _x.reset(new Eigen::VectorXd(_dataPoints.cols()));
    _y.reset(new Eigen::VectorXd(_dataPoints.rows()));

    //fill in the headers
    double x = _parameterBounds.getMinX();
    double * xIter = _x->data();
    for(int i = 0; i < nCol; ++i, ++xIter) {
        *xIter = x;
        x += resolution.getX();
    }
    *xIter = _parameterBounds.getMaxX();

    double y = _parameterBounds.getMinY();
    double * yIter = _y->data();
    for(int i = 0; i < nRow; ++i, ++yIter) {
        *yIter = y;
        y += resolution.getY();
    }
    *yIter = _parameterBounds.getMaxY();

    //compute the data points
    //outer loop over rows
    yIter = _y->data();
    for(int i = 0; i < _dataPoints.rows(); ++i, ++yIter) {
        xIter = _x->data();
        //inner loop over columns        
        for (int j = 0; j < _dataPoints.cols(); ++j, ++xIter) {
            //set grid point
            _dataPoints(i,j) = (*fillFunction)(*xIter, *yIter);
        }
    }
}

multifit::Cache::Functor::ConstPtr multifit::Cache::getRowFunctor(
    double const & y
) const {
    if (y < _parameterBounds.getMinY() || y >= _parameterBounds.getMaxY()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                y % _parameterBounds.getMinY() % _parameterBounds.getMaxY()).str()
        );
    }

    int i = static_cast<int>((y - _parameterBounds.getMinY()) / _yStep);
    double s = (y - (*_y)[i])/_getStep(i, _y);

    Eigen::VectorXd r(_dataPoints.row(i)*(1 - s) + _dataPoints.row(i+1)*s);
    return boost::make_shared<Cache::Functor const>(_x, r);
}

multifit::Cache::Functor::ConstPtr multifit::Cache::getRowDerivativeFunctor(
    double const & y
) const {
    if (y < _parameterBounds.getMinY() || y >= _parameterBounds.getMaxY()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                y % _parameterBounds.getMinY() % _parameterBounds.getMaxY()).str()
        );
    }

    int i = static_cast<int>((y - _parameterBounds.getMinY()) / _yStep);
    Eigen::VectorXd r((_dataPoints.row(i+1) - _dataPoints.row(i))/_getStep(i, _y));
    return boost::make_shared<Cache::Functor const>(_x, r);
}

multifit::Cache::Functor::ConstPtr multifit::Cache::getColFunctor(
    double const & x
) const {
    if (x < _parameterBounds.getMinX() || x >= _parameterBounds.getMaxX()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                x % _parameterBounds.getMinX() % _parameterBounds.getMaxX()).str()
        );
    }

    int i = static_cast<int>((x - _parameterBounds.getMinX()) / _xStep);
    double s = (x - (*_x)[i])/_getStep(i, _x);

    Eigen::VectorXd r(_dataPoints.col(i)*(1 - s) + _dataPoints.col(i+1)*s);
    return boost::make_shared<Cache::Functor const>(_y, r);
}

multifit::Cache::Functor::ConstPtr multifit::Cache::getColDerivativeFunctor(
    double const & x
) const {
    if (x < _parameterBounds.getMinX() || x >= _parameterBounds.getMaxX()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                x % _parameterBounds.getMinX() % _parameterBounds.getMaxX()).str()
        );
    }

    int i = static_cast<int>((x - _parameterBounds.getMinX()) / _xStep);
    Eigen::VectorXd r((_dataPoints.col(i+1) - _dataPoints.col(i))/_getStep(i, _x));
    return boost::make_shared<Cache::Functor const>(_y, r);
}

double multifit::Cache::Functor::operator() (double x) const {
    int i = static_cast<int>((x - (*_x)[0]) / _step);    
    if (i < 0 || i >= _x->size() - 1) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                x % (*_x)[0] % (*_x)[_x->size()-1]).str()
        );
    }

    double w1 = (x - (*_x)[i])/_getStep(i, _x);
    double w0 = 1 - w1;
    double p0 = _params[i];
    double p1 = _params[i+1];

    return w0*p0 + w1*p1;
}

double multifit::Cache::Functor::dParams(double const & x) const {
    int i = static_cast<int>((x - (*_x)[0]) / _step);    
    if (i < 0 || i >= _x->size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                x % (*_x)[0] % (*_x)[_x->size()-1]).str()
        );
    }
    double p0 = _params[i];
    double p1 = _params[i+1];
     
    return (p1 - p0)/ _getStep(i, _x);
}

multifit::Cache::Functor::Functor(
    boost::shared_ptr<Eigen::VectorXd> const & x, 
    std::vector<double> const & y
) : Base(y), _x(x) {
    initialize();  
}

multifit::Cache::Functor::Functor(
    boost::shared_ptr<Eigen::VectorXd> const & x,
    Eigen::VectorXd const & y 
) : Base(y.size()), _x(x) {
    initialize();
    std::copy(y.data(), y.data() + y.size(), _params.begin());
}

multifit::Cache::Functor::Base::Ptr multifit::Cache::Functor::clone() const {
    return boost::make_shared<Cache::Functor>(_x, _params);
}
void multifit::Cache::Functor::initialize() {
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

