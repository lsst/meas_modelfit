#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/afw/geom/ellipses.h"
namespace afwGeom = lsst::afw::geom;
namespace ellipses = lsst::afw::geom::ellipses;
namespace multifit = lsst::meas::multifit;

namespace Param {
    enum {FLUX, X, Y, GAMMA1=3,GAMMA2, KAPPA, N};
}
multifit::SmallGalaxyModelShape::SmallGalaxyModelShape(
        Eigen::VectorXd const & parameters,
        Eigen::MatrixXd const & covariance
) {
    set<X>(parameters[Param::X]);
    set<Y>(parameters[Param::Y]);

    double xVar = covariance(Param::X, Param::X);
    double yVar = covariance(Param::Y, Param::Y);
    set<X_ERR>(xVar*xVar);
    set<Y_ERR>(yVar*yVar);
    set<N>(parameters[Param::N]);

    ellipses::LogShear ls(
        parameters[Param::GAMMA1], 
        parameters[Param::GAMMA2], 
        parameters[Param::KAPPA]
    );

    //reparametrize portion of the covariance matrix in order to compute errors on
    //moments
    ellipses::Quadrupole q;
    ellipses::BaseCore::Jacobian mj = q.dAssign(ls);

    set<IXX>(q[ellipses::Quadrupole::IXX]);
    set<IXY>(q[ellipses::Quadrupole::IXY]);
    set<IYY>(q[ellipses::Quadrupole::IYY]);

    Eigen::Matrix<double, 3,3> momentsCov = mj.transpose() * 
            covariance.block<3,3>(Param::GAMMA1, Param::GAMMA1) * mj;       
    double ixxVar = momentsCov(ellipses::Quadrupole::IXX, ellipses::Quadrupole::IXX);
    double ixyVar = momentsCov(ellipses::Quadrupole::IXY, ellipses::Quadrupole::IXY);
    double iyyVar = momentsCov(ellipses::Quadrupole::IYY, ellipses::Quadrupole::IYY);

    set<IXX_ERR>(ixxVar*ixxVar);
    set<IXY_ERR>(ixyVar*ixyVar);
    set<IYY_ERR>(iyyVar*iyyVar); 
    
    //re-parametrize entrire covariance matrix for storage in e1,e2,r form
    ellipses::Distortion d;
    Eigen::Matrix<double, 7, 7> dj = Eigen::Matrix<double, 7, 7>::Identity();
    dj.block<3,3>(Param::GAMMA1, Param::GAMMA1) = d.dAssign(ls);
    Eigen::Matrix<double, 7, 7> cov = dj.transpose()*covariance*dj;
    
    set<E1>(d[ellipses::Distortion::E1]);
    set<E2>(d[ellipses::Distortion::E2]);
    set<R>(d[ellipses::Distortion::R]);
   
    set<FLUX_X_COV>(cov(Param::FLUX, Param::X));
    set<FLUX_Y_COV>(cov(Param::FLUX, Param::Y));
    set<FLUX_E1_COV>(cov(Param::FLUX, Param::GAMMA1));
    set<FLUX_E2_COV>(cov(Param::FLUX, Param::GAMMA2));
    set<FLUX_R_COV>(cov(Param::FLUX, Param::KAPPA));
    set<FLUX_N_COV>(cov(Param::FLUX, Param::N));

    set<X_X_COV>(cov(Param::X, Param::X));
    set<X_Y_COV>(cov(Param::X, Param::X));
    set<X_E1_COV>(cov(Param::X, Param::GAMMA1));
    set<X_E2_COV>(cov(Param::X, Param::GAMMA2));
    set<X_R_COV>(cov(Param::X, Param::KAPPA));
    set<X_N_COV>(cov(Param::X, Param::N));

    set<Y_Y_COV>(cov(Param::Y, Param::Y));
    set<Y_E1_COV>(cov(Param::Y, Param::GAMMA1));
    set<Y_E2_COV>(cov(Param::Y, Param::GAMMA2));
    set<Y_R_COV>(cov(Param::Y, Param::KAPPA));
    set<Y_N_COV>(cov(Param::Y, Param::N));

    set<E1_E1_COV>(cov(Param::GAMMA1, Param::GAMMA1));
    set<E1_E2_COV>(cov(Param::GAMMA1, Param::GAMMA2));
    set<E1_R_COV>(cov(Param::GAMMA1, Param::KAPPA));
    set<E1_N_COV>(cov(Param::GAMMA1, Param::N));

    set<E2_E2_COV>(cov(Param::GAMMA2, Param::GAMMA2));
    set<E2_R_COV>(cov(Param::GAMMA2, Param::KAPPA));
    set<E2_N_COV>(cov(Param::GAMMA2, Param::N));

    set<R_R_COV>(cov(Param::KAPPA, Param::KAPPA));
    set<R_N_COV>(cov(Param::KAPPA, Param::N));

    set<R_N_COV>(cov(Param::N, Param::N));
}

void multifit::SmallGalaxyModelShape::defineSchema(
    Schema::Ptr schema
) {
    Base::defineSchema(schema);
    schema->add(SchemaEntry("flux", FLUX, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e1", E1, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e2", E2, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("r", R, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("sersic", N, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxFluxCov", FLUX_FLUX_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxXCov", FLUX_X_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxYCov", FLUX_Y_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxE1Cov", FLUX_E1_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxE2Cov", FLUX_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxRCov", FLUX_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("fluxNCov", FLUX_N_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("xXCov", X_X_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("xYCov", X_Y_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("xE1Cov", X_E1_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("xE2Cov", X_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("xRCov", X_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("xNCov", X_N_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("yYCov", Y_Y_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("yE1Cov", Y_E1_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("yE2Cov", Y_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("yRCov", Y_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("yNCov", Y_N_COV, Schema::DOUBLE, 1));
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
