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

class DefaultInterpolationFunctionFactory : public multifit::InterpolationFunctionFactory {
public:

    static Registration registration;
    
    DefaultInterpolationFunctionFactory() : multifit::InterpolationFunctionFactory("", 0) {}

protected:

    virtual multifit::InterpolationFunction::ConstPtr execute(
        double const min, double const max, double const step,
        Eigen::VectorXd const & values
    ) const {
        return boost::make_shared<multifit::InterpolationFunction>(min, max, step, values);
    }
};

DefaultInterpolationFunctionFactory::Registration DefaultInterpolationFunctionFactory::registration(
    boost::make_shared<DefaultInterpolationFunctionFactory>()
);

} // unnamed namespace

multifit::InterpolationFunctionFactory::Registration::Registration(
    multifit::InterpolationFunctionFactory::Ptr const & p
) {
    getRegistry()[p->getName()] = p;
}

multifit::InterpolationFunctionFactory::ConstPtr multifit::InterpolationFunctionFactory::get(
    std::string const & name
) {
    Registry & registry = getRegistry();
    Registry::iterator i = registry.find(name);
    if (i == registry.end()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Interpolator not found in registry."
        );
    }
    return i->second;
}

multifit::InterpolationFunctionFactory::Registry & multifit::InterpolationFunctionFactory::getRegistry() {
    static Registry registry;
    return registry;
}

multifit::InterpolationFunctionFactory::InterpolationFunctionFactory(
    std::string const & name, int extraParameterSize
) : _name(name), _extraParameterSize(extraParameterSize) {}

multifit::InterpolationFunction::ConstPtr multifit::InterpolationFunctionFactory::operator()(
    double const min, double const max, double const step,
    Eigen::VectorXd const & values
) const {
    int n = int((max-min)/step) + 1;
    if (values.size() != n + getExtraParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            boost::str(
                boost::format(" Number of bins in fast dimension (%d) plus extra parameter size (%d) must match"
                              " length of vector values (%d).") % n % getExtraParameterSize() % values.size()
            )
        );
    }
    return execute(min, max, step, values);
}

std::map<std::string, multifit::Cache::Ptr> multifit::Cache::_registry;

multifit::Cache::ConstPtr multifit::Cache::get(std::string const & name) {
    return _registry[name];
}

multifit::Cache::ConstPtr multifit::Cache::make(
    double const slowMin, double const slowMax,
    double const fastMin, double const fastMax,
    int const nSlow, int const nFast,
    FillFunction const * fillFunction,
    InterpolatorFactory::ConstPtr const & interpolatorFactory,
    std::string const & name,
    bool const & doOverwrite
) {

    if(!doOverwrite && !name.empty()) {
        std::map<std::string, Ptr>::const_iterator i = _registry.find(name);
        if(i != _registry.end()){
            return i->second;
        }
    }
    Ptr cache(
        new Cache(
            slowMin, slowMax, fastMin, fastMax, 
            nSlow, nFast, 
            fillFunction, interpolatorFactory
        )
    );
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
    double const slowMin, double const slowMax,
    double const fastMin, double const fastMax,
    int const nSlow, int const nFast,
    FillFunction const * fillFunction,
    InterpolatorFactory::ConstPtr const & interpolatorFactory
) : _slowMin(slowMin), _slowMax(slowMax), _fastMin(fastMin), _fastMax(fastMax),
    _slowStep((slowMax-slowMin)/nSlow), _fastStep((fastMax-fastMin)/nFast),
    _interpolatorFactory(interpolatorFactory),
    _dataPoints(nSlow+1, nFast+1+interpolatorFactory->getExtraParameterSize())
{
    pexLog::Debug debug("lsst.meas.multifit.Cache");
    if(!interpolatorFactory) {
        _interpolatorFactory->get("");
    }

    //compute the data points
    //outer loop over slow dimension
    int nRows = nSlow+1, nCols = nFast+1;

    double slow=slowMin;
    for(int i = 0; i < nRows; ++i, slow+=_slowStep) {
        double fast = fastMin;
        debug.debug<5>("Filling row %d of %d.", i, nRows);
        //inner loop over fast dimension 
        for (int j = 0; j < nCols; ++j, fast += _fastStep) {
            //set grid point
            debug.debug<8>("Filling item %d of %d (column %d of %d).",
                           i*nCols + j, nCols * nRows, j, nCols);
            _dataPoints(i,j) = (*fillFunction)(slow, fast);
        }
    }

    slow = slowMin;
    for(int i = 0; i < nRows; ++i, slow += _slowStep) {
        debug.debug<5>("Filling extra parameters for row %d of %d.", i, nRows);
        _interpolatorFactory->fillExtraParameters(fastMax, slow, _dataPoints.row(i), *fillFunction);
    }
}
void multifit::Cache::save(std::string const & filepath) const {
    std::ofstream ofs(filepath.c_str());
    boost::archive::text_oarchive oa(ofs);
    oa << *this;
}

multifit::Cache::Interpolator::ConstPtr multifit::Cache::getInterpolator(
    double const slow
) const {
    if (slow < _slowMin || slow >= _slowMax) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                slow % _slowMin % _slowMax).str()
        );
    }

    int i = static_cast<int>((slow - _slowMin) / _slowStep);
    double s = (slow - (_slowMin+_slowStep*i))/_slowStep;

    Eigen::VectorXd r(_dataPoints.row(i)*(1 - s) + _dataPoints.row(i+1)*s);
    return (*_interpolatorFactory)(_fastMin, _fastMax, _fastStep, r);
}

multifit::Cache::Interpolator::ConstPtr multifit::Cache::getDerivativeInterpolator(
    double const slow
) const {
    if (slow < _slowMin || slow >= _slowMax) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                slow % _slowMin % _slowMax).str()
        );
    }

    int i = static_cast<int>((slow - _slowMin) / _slowStep);
    Eigen::VectorXd r((_dataPoints.row(i+1) - _dataPoints.row(i))/ _slowStep);
    return (*_interpolatorFactory)(_fastMin, _fastMax, _fastStep, r);
}

template <class Archive>
void multifit::Cache::serialize(Archive & ar, unsigned int const) {
    ar & _slowStep;
    ar & _fastStep;
    ar & _slowMin;
    ar & _slowMax;
    ar & _fastMin;
    ar & _fastMax;

    if (Archive::is_loading::value) {
        std::string interpolator;
        ar & interpolator;

        _interpolatorFactory = InterpolationFunctionFactory::get(interpolator);

        int nRow = getNSlow() + 1;
        int nCol = getNFast() + 1 + _interpolatorFactory->getExtraParameterSize();

        _dataPoints = Eigen::MatrixXd(nRow, nCol);
    }
    else {
        std::string interpolator = _interpolatorFactory->getName();
        ar & interpolator;
    }
    
    double * iData = _dataPoints.data();
    for(int i =0; i < _dataPoints.size(); ++i, ++iData) {
        ar & *iData;
    }
}

template void multifit::Cache::serialize<boost::archive::text_oarchive> (
    boost::archive::text_oarchive &, unsigned int const
);
template void multifit::Cache::serialize<boost::archive::text_iarchive> (
    boost::archive::text_iarchive &, unsigned int const
);

double multifit::InterpolationFunction::operator() (double point) const {
    if (point < _min || point >= _max) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                point % _min % _max).str()
        );
    }
    int i = static_cast<int>((point - _min) / _step);    


    double w1 = (point - (_min + _step*i))/_step;
    double w0 = 1 - w1;
    double v0 = _values[i];
    double v1 = _values[i+1];

    return w0*v0 + w1*v1;
}

double multifit::InterpolationFunction::d(double point) const {
    if (point < _min || point >= _max) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                point % _min % _max).str()
        );
    }
    int i = static_cast<int>((point - _min) / _step);    

    double v0 = _values[i];
    double v1 = _values[i+1];
     
    return (v1 - v0)/_step;
}
