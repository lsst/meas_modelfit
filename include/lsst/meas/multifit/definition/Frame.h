// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

class Frame {
public:
    
    Frame(
        ID id_, 
        Footprint::Ptr const & footprint_,
        lsst::ndarray::Array<Pixel,1,1> const & data_,
        lsst::ndarray::Array<Pixel,1,1> const & weights_ = lsst::ndarray::Array<Pixel,1,1>()
    ) : id(id_), filterId(Filter::UNKNOWN), 
        //calib(new Calib()), 
        wcs(), psf(),
        footprint(footprint_), data(data_), weights(weights_)
    {}

    Frame(
        ID id_,
        FilterId filterId_,
        //Calib::Ptr const & calib_,
        Wcs::Ptr const & wcs_,
        Psf::Ptr const & psf_,
        Footprint::Ptr const & footprint_,
        lsst::ndarray::Array<Pixel,1,1> const & data_,
        lsst::ndarray::Array<Pixel,1,1> const & weights_
    ) : id(id_), filterId(filterId_), 
        //calib(calib_), 
        wcs(wcs_), psf(psf_), footprint(footprint_), 
        data(data_), weights(weights_)
    {}

    Frame(Frame const & other) :
        id(other.id), filterId(other.filterId), 
        //calib(other.calib),
        wcs(other.wcs), psf(other.psf), footprint(other.footprint), 
        data(other.data), weights(other.weights)
    {}

    ID const id;
    
    FilterId filterId;
    //Calib::Ptr calib;    
    Wcs::Ptr wcs;
    Psf::Ptr psf;

    Footprint::Ptr footprint;
    lsst::ndarray::Array<Pixel,1,1> data;
    lsst::ndarray::Array<Pixel,1,1> weights;
};

std::ostream & operator<<(std::ostream & os, Frame const & frame);

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Frame
