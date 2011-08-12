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

#ifndef LSST_MEAS_MULTIFIT_Definition
#define LSST_MEAS_MULTIFIT_Definition

#include <fstream>
#include "lsst/meas/multifit/definition/Frame.h"
#include "lsst/meas/multifit/definition/ObjectComponent.h"
#include "lsst/meas/multifit/containers/Set.h"
#include "lsst/afw/formatters/WcsFormatter.h"
#include "boost/serialization/nvp.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"

namespace boost {
namespace serialization {
    class access;
}}

namespace lsst { namespace meas { namespace multifit { namespace definition {



class Definition {
public:    
    typedef definition::Frame Frame;
    typedef definition::ObjectComponent ObjectComponent;

    typedef containers::MutableSet<Frame> FrameSet;
    typedef containers::MutableSet<ObjectComponent> ObjectComponentSet; 

    Definition() {}

    explicit Definition(Wcs::Ptr const & wcs) : _wcs(wcs) {}

    Definition(Definition const & other);

    FrameSet frames;
    ObjectComponentSet objects;

    Wcs::Ptr getWcs() const { return _wcs; }
    static Definition load(std::string const & file) {
        Definition definition;
        std::ifstream ifs(file.c_str());
        boost::archive::text_iarchive ar(ifs);
        ar & definition;
        return definition;
    }
    
    void save(std::string const & file){
        std::ofstream ofs(file.c_str());
        boost::archive::text_oarchive ar(ofs);

        ar& *this;
    };

private:
    Wcs::Ptr _wcs;

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, unsigned int const version) {
        ar & boost::serialization::make_nvp("frames", frames);
        ar & boost::serialization::make_nvp("objects", objects);
        ar & boost::serialization::make_nvp("wcs", _wcs);
    }
};

} // namespace definition

typedef definition::Definition Definition;

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Definition
