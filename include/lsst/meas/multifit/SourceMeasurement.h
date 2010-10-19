#ifndef LSST_MEAS_MULTIFIT_SOURCE_MEASUREMENT_H
#define LSST_MEAS_MULTIFIT_SOURCE_MEASUREMENT_H

#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Shape.h"
#include <Eigen/Core>

namespace lsst {
namespace meas {
namespace multifit {

class PointSourceModelPhotometry : public lsst::afw::detection::Photometry {
public:
    typedef lsst::afw::detection::Photometry Base;
    virtual ~PointSourceModelPhotometry(){}

    PointSourceModelPhotometry(double flux, double fluxErr) : Base(flux, fluxErr) {}
#ifndef SWIG
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema) {
        Base::defineSchema(schema);
        schema->setComponent("psModel");
    }
#endif
};

class SmallGalaxyModelPhotometry : public lsst::afw::detection::Photometry {
public:
    typedef lsst::afw::detection::Photometry Base;

    virtual ~SmallGalaxyModelPhotometry(){}
    SmallGalaxyModelPhotometry(double flux, double fluxErr) : Base(flux, fluxErr) {}
#ifndef SWIG
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema) {
        Base::defineSchema(schema);
        schema->setComponent("sgModel");
    }
#endif
};

class SmallGalaxyModelShape : public lsst::afw::detection::Shape {
public:
    typedef lsst::afw::detection::Schema Schema;
    typedef lsst::afw::detection::SchemaEntry SchemaEntry;
    typedef lsst::afw::detection::Shape Base;
    typedef lsst::afw::detection::Measurement<Base> Measurement;

    enum {N=Base::NVALUE, FLUX, E1, E2, R,
        FLUX_FLUX_COV, FLUX_X_COV, FLUX_Y_COV, FLUX_E1_COV, FLUX_E2_COV, FLUX_R_COV, FLUX_N_COV, 
        X_X_COV, X_Y_COV, X_E1_COV, X_E2_COV, X_R_COV, X_N_COV,
        Y_Y_COV, Y_E1_COV, Y_E2_COV, Y_R_COV, Y_N_COV,
        E1_E1_COV, E1_E2_COV, E1_R_COV, E1_N_COV,
        E2_E2_COV, E2_R_COV, E2_N_COV,
        R_R_COV, R_N_COV,
        N_N_COV, NVALUE
    };

    SmallGalaxyModelShape(
        Eigen::VectorXd const & parameters,
        Eigen::MatrixXd const & covariance
    );
    virtual ~SmallGalaxyModelShape(){}
#ifndef SWIG
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema);
#endif
};


}}}

#endif
