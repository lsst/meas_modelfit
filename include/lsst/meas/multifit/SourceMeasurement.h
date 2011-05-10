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

private:
    PointSourceModelPhotometry(void) : lsst::afw::detection::Photometry() { }
    LSST_SERIALIZE_PARENT(lsst::afw::detection::Photometry);
};

class SmallGalaxyModelPhotometry : public lsst::afw::detection::Photometry {
public:
    typedef lsst::afw::detection::Schema Schema;
    typedef lsst::afw::detection::SchemaEntry SchemaEntry;
    typedef lsst::afw::detection::Photometry Base;
    typedef lsst::afw::detection::Measurement<Base> Measurement;

    enum {AMPLITUDE=Base::NVALUE, E1, E2, R, 
        AMP_AMP_COV, AMP_E1_COV, AMP_E2_COV, AMP_R_COV, 
        E1_E1_COV, E1_E2_COV, E1_R_COV,  
        E2_E2_COV, E2_R_COV,
        R_R_COV, 
        NVALUE
    };

    SmallGalaxyModelPhotometry(
        Eigen::VectorXd const & parameters,
        Eigen::MatrixXd const & covariance
    );
    SmallGalaxyModelPhotometry(
        std::vector<double> const & parameters,
        Eigen::MatrixXd const & covariance
    );

    virtual ~SmallGalaxyModelPhotometry(){}
#ifndef SWIG
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema) {
        Base::defineSchema(schema);
        schema->add(SchemaEntry("amplitude", AMPLITUDE, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e1", E1, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e2", E2, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("r", R, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("ampAmpCov", AMP_AMP_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("ampE1Cov", AMP_E1_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("ampE2Cov", AMP_E2_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("ampRCov", AMP_R_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e1E1Cov", E1_E1_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e1E2Cov", E1_E2_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e1RCov", E1_R_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e2E2Cov", E2_E2_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("e2RCov", E2_R_COV, Schema::DOUBLE, 1));
        schema->add(SchemaEntry("rRCov", R_R_COV, Schema::DOUBLE, 1));
        schema->setComponent("sgModel"); 
    }
#endif

private:
    void fill(
        double const * parameters,
        Eigen::MatrixXd const & covariance
    );
    SmallGalaxyModelPhotometry(void) : lsst::afw::detection::Photometry() { }
    LSST_SERIALIZE_PARENT(lsst::afw::detection::Photometry);
};

}}}

LSST_REGISTER_SERIALIZER(lsst::meas::multifit::PointSourceModelPhotometry);
LSST_REGISTER_SERIALIZER(lsst::meas::multifit::SmallGalaxyModelPhotometry);

#endif
