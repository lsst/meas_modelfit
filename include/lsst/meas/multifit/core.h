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
 
/**
 * @file
 * Collection of useful typedefs used throughout the lsst::meas::multifit
 * namespace
 */
#ifndef LSST_MEAS_MULTIFIT_CORE_H
#define LSST_MEAS_MULTIFIT_CORE_H

#include <Eigen/Core>

#include "lsst/afw/geom.h"
#include "lsst/afw/geom/deprecated.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace lsst {
namespace meas {
namespace multifit {

typedef double Pixel;
typedef double Parameter;

typedef Parameter * ParameterIterator;
typedef Parameter const * ParameterConstIterator;
typedef std::map<int, Parameter> ParameterMap;

typedef Eigen::aligned_allocator<char> Allocator;
typedef Eigen::VectorXd ParameterVector;
typedef Eigen::Map<Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> > MatrixMap;
typedef Eigen::Block<MatrixMap> MatrixMapBlock;
typedef Eigen::Map<Eigen::Matrix<Pixel, Eigen::Dynamic, 1> > VectorMap;


class ParameterRangeException: public lsst::pex::exceptions::RangeErrorException {
public:
    ParameterRangeException(
        char const * file, int line, char const * func,
        std::string const & message
    ) : lsst::pex::exceptions::RangeErrorException(file, line, func, message) {}
    ParameterRangeException(
        char const * file, int line, char const * func,
        std::string const & message,
        ParameterMap outOfRangeParameterMap
    ) : lsst::pex::exceptions::RangeErrorException(file, line, func, message),
        _map(outOfRangeParameterMap) 
    {}
    
    virtual char const * getType(void) const throw() {
        return "lsst::meas::multifit::ParameterRangeException *";
    }
    virtual lsst::pex::exceptions::Exception* clone(void) const {
        return new ParameterRangeException(*this);
    }
    ParameterMap const & getOutOfRangeParameterMap() const {
        return _map;
    }
    bool notInRange(int const parameter) const {
        return _map.find(parameter) == _map.end();
    }
private:
    ParameterMap _map;
};

}}} // namespace lsst::meas::multifit

#endif
