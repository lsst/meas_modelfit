#include "lsst/meas/multifit/SourceMeasurement.h"

#include "lsst/afw/geom/ellipses.h"

namespace lsst { namespace meas { namespace multifit {

template <int nCoeff>
void ShapeletModelPhotometry<nCoeff>::defineSchema(lsst::afw::detection::Schema::Ptr schema) {
    schema->clear();
    schema->add(afw::detection::SchemaEntry("flux",    FLUX,     afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("fluxErr", FLUX_ERR, afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("e1",      E1,       afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("e2",      E2,       afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("radius",  RADIUS,   afw::detection::Schema::DOUBLE, 1,
                                            "pixels"));
    schema->add(afw::detection::SchemaEntry("coefficients", COEFFICIENTS, afw::detection::Schema::DOUBLE,
                                            nCoeff));
}

}}}
