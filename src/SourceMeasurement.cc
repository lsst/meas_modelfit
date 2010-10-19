#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/afw/geom/ellipses.h"
namespace afwGeom = lsst::afw::geom;
namespace ellipses = lsst::afw::geom::ellipses;
namespace multifit = lsst::meas::multifit;

namespace Param {
    enum {FLUX=0, GAMMA1,GAMMA2, KAPPA, N, NPARAM};
}
multifit::SmallGalaxyModelPhotometry::SmallGalaxyModelPhotometry(
        std::vector<double> const & parameters,
        Eigen::MatrixXd const & covariance
) {
    ellipses::LogShear ls(
        parameters[Param::GAMMA1], 
        parameters[Param::GAMMA2], 
        parameters[Param::KAPPA]
    );

    //re-parametrize entrire covariance matrix for storage in e1,e2,r form
    ellipses::Distortion d;
    Eigen::Matrix<double, Param::NPARAM, Param::NPARAM> dj = 
        Eigen::Matrix<double, Param::NPARAM, Param::NPARAM>::Identity();
    dj.block<3,3>(Param::GAMMA1, Param::GAMMA1) = d.dAssign(ls);
    Eigen::Matrix<double, Param::NPARAM, Param::NPARAM> cov = 
        dj.transpose()*covariance*dj;
   
    set<FLUX>(parameters[Param::FLUX]);
    set<E1>(d[ellipses::Distortion::E1]);
    set<E2>(d[ellipses::Distortion::E2]);
    set<R>(d[ellipses::Distortion::R]);
    set<N>(parameters[Param::N]);
  
    double fluxFluxCov = cov(Param::FLUX, Param::FLUX);
    set<FLUX_FLUX_COV>(fluxFluxCov);
    set<FLUX_ERR>(sqrt(fluxFluxCov));
    set<FLUX_E1_COV>(cov(Param::FLUX, Param::GAMMA1));
    set<FLUX_E2_COV>(cov(Param::FLUX, Param::GAMMA2));
    set<FLUX_R_COV>(cov(Param::FLUX, Param::KAPPA));
    set<FLUX_N_COV>(cov(Param::FLUX, Param::N));

    set<E1_E1_COV>(cov(Param::GAMMA1, Param::GAMMA1));
    set<E1_E2_COV>(cov(Param::GAMMA1, Param::GAMMA2));
    set<E1_R_COV>(cov(Param::GAMMA1, Param::KAPPA));
    set<E1_N_COV>(cov(Param::GAMMA1, Param::N));

    set<E2_E2_COV>(cov(Param::GAMMA2, Param::GAMMA2));
    set<E2_R_COV>(cov(Param::GAMMA2, Param::KAPPA));
    set<E2_N_COV>(cov(Param::GAMMA2, Param::N));

    set<R_R_COV>(cov(Param::KAPPA, Param::KAPPA));
    set<R_N_COV>(cov(Param::KAPPA, Param::N));

    set<N_N_COV>(cov(Param::N, Param::N));
}

void multifit::SmallGalaxyModelPhotometry::defineSchema(
    Schema::Ptr schema
) {
    Base::defineSchema(schema);
    schema->add(SchemaEntry("e1", E1, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e2", E2, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("r", R, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("sersic", N, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxFluxCov", FLUX_FLUX_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxE1Cov", FLUX_E1_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxE2Cov", FLUX_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxRCov", FLUX_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxNCov", FLUX_N_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e1E1Cov", E1_E1_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e1E2Cov", E1_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e1RCov", E1_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e1NCov", E1_N_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e2E2Cov", E2_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e2RCov", E2_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e2NCov", E2_N_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("rRCov", R_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("rNCov", R_N_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("nNCov", N_N_COV, Schema::DOUBLE, 1));

    schema->setComponent("sgModel");
}
