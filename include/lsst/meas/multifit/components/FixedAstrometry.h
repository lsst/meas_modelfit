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
 
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "lsst/afw/image/Utils.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/components/Astrometry.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class FixedAstrometry : public Astrometry {
public:
    typedef boost::shared_ptr<FixedAstrometry> Ptr;
    typedef boost::shared_ptr<FixedAstrometry const> ConstPtr;
    
    explicit FixedAstrometry(lsst::afw::geom::Point2D const & position) 
      : Astrometry(position) {}

    explicit FixedAstrometry(lsst::afw::coord::Coord const & coord) 
      : Astrometry(coord) {}
    virtual ~FixedAstrometry() {}
    
    virtual DerivativeMatrix const & differentiate() const {
        static DerivativeMatrix i;
        return i;
    }
protected:
    virtual int const getParameterSize() const { return 0; }

    virtual Astrometry::Ptr create(
        boost::shared_ptr<ParameterVector const> const &,
        size_t const &
    ) const {
        return boost::make_shared<FixedAstrometry>(computePosition());
    }
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H
