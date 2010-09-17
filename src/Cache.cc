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
 
#include "lsst/meas/multifit/Cache.h"

#include <cmath>
#include <boost/make_shared.hpp>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <fstream>
#include <boost/serialization/list.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "lsst/pex/exceptions/Runtime.h"

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"

namespace multifit = lsst::meas::multifit;
namespace pexLog = lsst::pex::logging;

namespace {
    //helper function to determine the difference between element i+1 and i
    //of vetor r
    double const getStep(int const &i, Eigen::VectorXd const& r) {
        return r[i+1] - r[i];
    }
}
std::map<std::string, multifit::Cache::Ptr> multifit::Cache::_registry;

multifit::Cache::ConstPtr multifit::Cache::get(std::string const & name) {
    return _registry[name];
}

multifit::Cache::ConstPtr multifit::Cache::make(
    lsst::afw::geom::BoxD const & parameterBounds,
    lsst::afw::geom::Extent2D const & resolution,
    FillFunction const * fillFunction,
    std::string const & name,
    bool const & doOverwrite
) {

    if(!doOverwrite && !name.empty()) {
        std::map<std::string, Ptr>::const_iterator i = _registry.find(name);
        if(i != _registry.end()){
            return i->second;
        }
    }
    Ptr cache(new Cache(parameterBounds, resolution, fillFunction));

    if(doOverwrite && !name.empty()) {
            _registry[name] = cache;
    }
    
    return cache;
}

multifit::Cache::ConstPtr multifit::Cache::load(
    std::string const & filepath,
    std::string const & name,
    bool const & doOverwrite
) {

    if(!doOverwrite && !name.empty()) {
        std::map<std::string, Ptr>::const_iterator i = _registry.find(name);
        if(i != _registry.end()){
            return i->second;
        }
    }
    Ptr cache(new Cache());
    std::ifstream ifs(filepath.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> *cache;
    
    if(doOverwrite && !name.empty()) {
            _registry[name] = cache;
    }

    return cache;
}


multifit::Cache::Cache(
    lsst::afw::geom::BoxD const & parameterBounds,
    lsst::afw::geom::Extent2D const & resolution,
    FillFunction const * fillFunction
) : _parameterBounds(parameterBounds),
    _xStep(resolution.getX()),
    _yStep(resolution.getY())    
{
    pexLog::Debug debug("lsst.meas.multifit.Cache");

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
    _x = Eigen::VectorXd(_dataPoints.cols());
    _y = Eigen::VectorXd(_dataPoints.rows());

    //fill in the headers
    double x = _parameterBounds.getMinX();
    double * xIter = _x.data();
    for(int i = 0; i < nCol; ++i, ++xIter) {
        *xIter = x;
        x += resolution.getX();
    }
    *xIter = _parameterBounds.getMaxX();

    double y = _parameterBounds.getMinY();
    double * yIter = _y.data();
    for(int i = 0; i < nRow; ++i, ++yIter) {
        *yIter = y;
        y += resolution.getY();
    }
    *yIter = _parameterBounds.getMaxY();

    //compute the data points
    //outer loop over rows
    yIter = _y.data();
    for(int i = 0; i < _dataPoints.rows(); ++i, ++yIter) {
        xIter = _x.data();
        debug.debug<5>("Filling row %d of %d.", i, _dataPoints.rows());
        //inner loop over columns        
        for (int j = 0; j < _dataPoints.cols(); ++j, ++xIter) {
            //set grid point
            debug.debug<8>("Filling item %d of %d.", i*_dataPoints.cols() + j, _dataPoints.size());
            _dataPoints(i,j) = (*fillFunction)(*xIter, *yIter);
        }
    }
}
void multifit::Cache::save(std::string const & filepath) const {
    std::ofstream ofs(filepath.c_str());
    boost::archive::text_oarchive oa(ofs);
    oa << *this;
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
    double s = (y - _y[i])/::getStep(i, _y);

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
    Eigen::VectorXd r((_dataPoints.row(i+1) - _dataPoints.row(i))/::getStep(i, _y));
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
    double s = (x - _x[i])/::getStep(i, _x);

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
    Eigen::VectorXd r((_dataPoints.col(i+1) - _dataPoints.col(i))/::getStep(i, _x));
    return boost::make_shared<Cache::Functor const>(_y, r);
}

template <class Archive>
void multifit::Cache::serialize(Archive & ar, unsigned int const) {
    ar & _xStep;
    ar & _yStep;

    if (Archive::is_loading::value) {
        double x, y, width, height;        
        ar & x;
        ar & y;
        ar & width;
        ar & height;

        lsst::afw::geom::Point2D min = lsst::afw::geom::Point2D::make(x, y);
        lsst::afw::geom::Extent2D dim = lsst::afw::geom::Extent2D::make(width, height);
        _parameterBounds = lsst::afw::geom::BoxD(min, dim);


        int nCol = width / _xStep+1;
        int nRow = height / _yStep+1;

        _dataPoints = Eigen::MatrixXd(nRow, nCol);
        _x = Eigen::VectorXd(nCol);
        _y = Eigen::VectorXd(nRow);
    }
    else {
        double x = _parameterBounds.getMinX();
        double y = _parameterBounds.getMinY();
        double width = _parameterBounds.getWidth();
        double height = _parameterBounds.getHeight();       

        ar & x;
        ar & y;
        ar & width;
        ar & height;
    }
    
    double * iData = _dataPoints.data();
    for(int i =0; i < _dataPoints.size(); ++i, ++iData) {
        ar & *iData;
    }
    double *iX = _x.data();
    for (int i =0; i < _x.size(); ++i, ++iX) {
        ar & *iX;
    }
    double *iY = _y.data();
    for (int i =0; i < _y.size(); ++i, ++iY) {
        ar & *iY;
    }
}

template void multifit::Cache::serialize<boost::archive::text_oarchive> (
    boost::archive::text_oarchive &, unsigned int const
);
template void multifit::Cache::serialize<boost::archive::text_iarchive> (
    boost::archive::text_iarchive &, unsigned int const
);
double multifit::InterpolationFunction::operator() (double x) const {
    int i = static_cast<int>((x - _x[0]) / _step);    
    if (i < 0 || i >= _x.size() - 1) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                x % _x[0] % _x[_x.size()-1]).str()
        );
    }

    double w1 = (x - _x[i])/::getStep(i, _x);
    double w0 = 1 - w1;
    double p0 = _params[i];
    double p1 = _params[i+1];

    return w0*p0 + w1*p1;
}

double multifit::InterpolationFunction::dParams(double const & x) const {
    int i = static_cast<int>((x - _x[0]) / _step);    
    if (i < 0 || i >= _x.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                x % _x[0] % _x[_x.size()-1]).str()
        );
    }
    double p0 = _params[i];
    double p1 = _params[i+1];
     
    return (p1 - p0)/ ::getStep(i, _x);
}

multifit::InterpolationFunction::InterpolationFunction(
    Eigen::VectorXd const & x, 
    std::vector<double> const & y
) : Base(y), _x(x) {
    initialize();  
}

multifit::InterpolationFunction::InterpolationFunction(
    Eigen::VectorXd const & x,
    Eigen::VectorXd const & y 
) : Base(y.size()), _x(x) {
    initialize();
    std::copy(y.data(), y.data() + y.size(), _params.begin());
}

multifit::InterpolationFunction::Base::Ptr multifit::Cache::Functor::clone() const {
    return boost::make_shared<InterpolationFunction>(_x, _params);
}
void multifit::InterpolationFunction::initialize() {
    if(_x.size() < 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Not enough points. At least 2 data points needed"
        );
    }
    if(getNParameters() != static_cast<unsigned int>(_x.size())) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Length of vector x must match length of vector y"
        );
    }
    _step = _x[1] - _x[0];
}

