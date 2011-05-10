#include "lsst/meas/multifit/SourceMeasurement.h"

#include "lsst/afw/geom/ellipses.h"

namespace afwGeom = lsst::afw::geom;
namespace ellipses = lsst::afw::geom::ellipses;
namespace multifit = lsst::meas::multifit;

namespace Param {
    enum {AMPLITUDE=0, GAMMA1,GAMMA2, KAPPA, NPARAM};
}

void multifit::SmallGalaxyModelPhotometry::fill(
        double const * parameters,
        Eigen::MatrixXd const & covariance
) {
    init();

    /*
    ellipses::LogShear ls(
        parameters[Param::GAMMA1], 
        parameters[Param::GAMMA2], 
        parameters[Param::KAPPA]
    );

    //re-parametrize entire covariance matrix for storage in e1,e2,r form
    ellipses::Distortion d;
    Eigen::Matrix<double, Param::NPARAM, Param::NPARAM> dj = 
        Eigen::Matrix<double, Param::NPARAM, Param::NPARAM>::Identity();
    dj.block<3,3>(Param::GAMMA1, Param::GAMMA1) = d.dAssign(ls);
    Eigen::Matrix<double, Param::NPARAM, Param::NPARAM> cov = 
        dj.transpose()*covariance*dj;
  
    double r = d[ellipses::Distortion::R];
    double amp = parameters[Param::AMPLITUDE];
    set<AMPLITUDE>(amp);
    set<E1>(d[ellipses::Distortion::E1]);
    set<E2>(d[ellipses::Distortion::E2]);
    set<R>(d[ellipses::Distortion::R]);


    set<FLUX>(amp*r*r);
  
    double ampVar = cov(Param::AMPLITUDE, Param::AMPLITUDE);
    double radiusVar = cov(Param::KAPPA, Param::KAPPA);
    set<AMP_AMP_COV>(ampVar);

    double dfDr = 2*p*amp*r;
    double dfDa = p*r*r;
    double fluxVar = radiusVar * dfDr*dfDr + ampVar*dfDa*dfDa;

    set<FLUX_ERR>(fluxVar);
    set<AMP_E1_COV>(cov(Param::AMPLITUDE, Param::GAMMA1));
    set<AMP_E2_COV>(cov(Param::AMPLITUDE, Param::GAMMA2));
    set<AMP_R_COV>(cov(Param::AMPLITUDE, Param::KAPPA));
    set<E1_E1_COV>(cov(Param::GAMMA1, Param::GAMMA1));
    set<E1_E2_COV>(cov(Param::GAMMA1, Param::GAMMA2));
    set<E1_R_COV>(cov(Param::GAMMA1, Param::KAPPA));
    set<E2_E2_COV>(cov(Param::GAMMA2, Param::GAMMA2));
    set<E2_R_COV>(cov(Param::GAMMA2, Param::KAPPA));
    set<R_R_COV>(radiusVar);
    */
}


