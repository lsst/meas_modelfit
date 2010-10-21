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

struct Flags{
    enum {
        FAIL_INIT_PS_NAN        = 0x00010000, 
        FAIL_INIT_PS_NO_PIX     = 0x00020000, 
        FAIL_FIT_PS             = 0x00040000, 
        FAIL_INIT_SG_NAN        = 0x00080000, 
        FAIL_INIT_SG_MOMENTS    = 0x00100000, 
        FAIL_INIT_SG_NO_PIX     = 0x00200000, 
        FAIL_FIT_SG_DISTORTION  = 0x00400000, 
        FAIL_FIT_SG_SERSIC      = 0x00800000, 
        FAIL_FIT_SG_OPTIMIZER   = 0x01000000,
        PS_MAX_ITERATIONS       = 0x02000000,
        SG_MAX_ITERATIONS       = 0x04000000,

        FAIL_INIT_PS = FAIL_INIT_PS_NAN | FAIL_INIT_PS_NO_PIX,
        FAIL_PS = FAIL_INIT_PS | FAIL_FIT_PS,
        FAIL_INIT_SG= FAIL_INIT_SG_NAN | FAIL_INIT_SG_MOMENTS | FAIL_INIT_SG_NO_PIX,
        FAIL_FIT_SG = FAIL_FIT_SG_DISTORTION|FAIL_FIT_SG_SERSIC|FAIL_FIT_SG_OPTIMIZER,
        FAIL_SG=FAIL_INIT_SG|FAIL_FIT_SG
    };
};
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
    typedef lsst::afw::detection::Schema Schema;
    typedef lsst::afw::detection::SchemaEntry SchemaEntry;
    typedef lsst::afw::detection::Photometry Base;
    typedef lsst::afw::detection::Measurement<Base> Measurement;

    enum {E1=Base::NVALUE, E2, R, N, 
        FLUX_FLUX_COV, FLUX_E1_COV, FLUX_E2_COV, FLUX_R_COV, FLUX_N_COV, 
        E1_E1_COV, E1_E2_COV, E1_R_COV, E1_N_COV,
        E2_E2_COV, E2_R_COV, E2_N_COV,
        R_R_COV, R_N_COV,
        N_N_COV, 
        NVALUE
    };

    SmallGalaxyModelPhotometry(
        std::vector<double> const & parameters,
        Eigen::MatrixXd const & covariance
    );
    virtual ~SmallGalaxyModelPhotometry(){}
#ifndef SWIG
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema);
#endif
};


}}}

#endif
