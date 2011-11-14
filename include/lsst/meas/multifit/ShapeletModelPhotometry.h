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

#ifndef LSST_MEAS_MULTIFIT_ShapeletModelPhotometry_H
#define LSST_MEAS_MULTIFIT_ShapeletModelPhotometry_H

#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/ndarray.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/detection/Peak.h"
#include "lsst/meas/multifit/ModelBasis.h"
#include "lsst/meas/multifit/SourceMeasurement.h"

namespace lsst {
namespace meas {
namespace multifit {

class ShapeletModelPhotometry : public lsst::afw::detection::Photometry {
public:
    typedef lsst::afw::detection::Schema Schema;
    typedef lsst::afw::detection::SchemaEntry SchemaEntry;
    typedef lsst::afw::detection::Photometry Base;
    typedef lsst::afw::detection::Measurement<Base> Measurement;

    enum {
        FLUX = Base::FLUX,
        FLUX_ERR,
        STATUS,
        E1, E2, RADIUS, 
        COEFFICIENTS
    };

    ShapeletModelPhotometry(boost::int64_t flag, int nCoeff=0);
    ShapeletModelPhotometry(boost::int64_t flag, int nCoeff,
        double flux, double fluxErr,
        double e1, double e2, double radius,
        ndarray::Array<double const, 1,1> coefficients
    );

    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema);
    virtual boost::int64_t getFlag() const {return get<STATUS, boost::int64_t>(); }

private:
    int _nCoeff;

    ShapeletModelPhotometry() : lsst::afw::detection::Photometry() { init(); }
    LSST_SERIALIZE_PARENT(lsst::afw::detection::Photometry);
};

}}}

LSST_REGISTER_SERIALIZER(lsst::meas::multifit::ShapeletModelPhotometry);

#endif
