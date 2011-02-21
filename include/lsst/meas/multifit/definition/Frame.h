// -*- LSST-C++ -*-
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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_Frame
#define LSST_MEAS_MULTIFIT_DEFINITION_Frame

#include "lsst/ndarray.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit { namespace definition {

class Filter {
public:

    typedef boost::shared_ptr<Filter> Ptr;

    static Filter * get(std::string const & name);

    bool operator==(Filter const & other) const { return this == &other; }
    bool operator!=(Filter const & other) const { return !(*this == other); }

    std::string const name;

private:
    Filter(std::string const & name_) : name(name_) {}
};

class Frame {
public:
    
    Frame(
        ID id_, 
        Footprint::Ptr const & footprint_,
        ndarray::Array<double,1,1> const & data_,
        ndarray::Array<double,1,1> const & weights_ = ndarray::Array<double,1,1>()
    ) : id(id_), filter(), startTime(0), stopTime(0), wcs(), psf(),
        footprint(footprint_), data(data_), weights(weights_)
    {}

    Frame(
        ID id_,
        Filter const * filter_,
        Timestamp startTime_,
        Timestamp stopTime_,
        Wcs::Ptr const & wcs_,
        Psf::Ptr const & psf_,
        Footprint::Ptr const & footprint_,
        ndarray::Array<double,1,1> const & data_,
        ndarray::Array<double,1,1> const & weights_
    ) : id(id_), filter(filter_), startTime(startTime_), stopTime(stopTime_), 
        wcs(wcs_), psf(psf_), footprint(footprint_), data(data_), weights(weights_)
    {}

    Frame(Frame const & other) :
        id(other.id), filter(other.filter), startTime(other.startTime), stopTime(other.stopTime),
        wcs(other.wcs), psf(other.psf), footprint(other.footprint), data(other.data), weights(other.weights)
    {}

    ID const id;
    
    Filter const * filter;

    Timestamp startTime;
    Timestamp stopTime;

    Wcs::Ptr wcs;
    Psf::Ptr psf;

    Footprint::Ptr footprint;
    ndarray::Array<double,1,1> data;
    ndarray::Array<double,1,1> weights;
};

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Frame
