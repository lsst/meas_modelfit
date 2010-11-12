#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/ModifiedSersic.h"
#include "lsst/afw/geom/ellipses.h"

namespace afwGeom = lsst::afw::geom;
namespace ellipses = lsst::afw::geom::ellipses;
namespace multifit = lsst::meas::multifit;

namespace Param {
    enum {AMPLITUDE=0, GAMMA1,GAMMA2, KAPPA, N, NPARAM};
}

multifit::SmallGalaxyModelPhotometry::SmallGalaxyModelPhotometry(
        Eigen::VectorXd const & parameters,
        Eigen::MatrixXd const & covariance,
        double innerSersicRadius, double outerSersicRadius
) {
    fill(parameters.data(), covariance, innerSersicRadius, outerSersicRadius);
}
multifit::SmallGalaxyModelPhotometry::SmallGalaxyModelPhotometry(
        std::vector<double> const & parameters,
        Eigen::MatrixXd const & covariance,
        double innerSersicRadius, double outerSersicRadius
) {
    fill(&parameters[0], covariance, innerSersicRadius, outerSersicRadius);
}

void multifit::SmallGalaxyModelPhotometry::fill(
        double const * parameters,
        Eigen::MatrixXd const & covariance,
        double innerSersicRadius, double outerSersicRadius
) {
    init();

    ellipses::LogShear ls(
        parameters[Param::GAMMA1], 
        parameters[Param::GAMMA2], 
        parameters[Param::KAPPA]
    );

    SersicCache::ConstPtr cache = components::SersicMorphology::getSersicCache();
    double sersicIndex = cache->convertParameterToSersic(parameters[Param::N]);

    //re-parametrize entire covariance matrix for storage in e1,e2,r form
    ellipses::Distortion d;
    Eigen::Matrix<double, Param::NPARAM, Param::NPARAM> dj = 
        Eigen::Matrix<double, Param::NPARAM, Param::NPARAM>::Identity();
    dj(Param::N, Param::N) = cache->differentiateParameterToSersic(parameters[Param::N]);
    dj.block<3,3>(Param::GAMMA1, Param::GAMMA1) = d.dAssign(ls);
    Eigen::Matrix<double, Param::NPARAM, Param::NPARAM> cov = 
        dj.transpose()*covariance*dj;
  
    double r = d[ellipses::Distortion::R];
    double amp = parameters[Param::AMPLITUDE];
    set<AMPLITUDE>(amp);
    set<E1>(d[ellipses::Distortion::E1]);
    set<E2>(d[ellipses::Distortion::E2]);
    set<R>(d[ellipses::Distortion::R]);

    set<N>(sersicIndex);
    ModifiedSersicFunction func(
        sersicIndex, innerSersicRadius, outerSersicRadius
    );
    double p=func.integrate(outerSersicRadius, 1)*2.0*M_PI;
    set<FLUX>(p*amp*r*r);
  
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
    set<AMP_N_COV>(cov(Param::AMPLITUDE, Param::N));

    set<E1_E1_COV>(cov(Param::GAMMA1, Param::GAMMA1));
    set<E1_E2_COV>(cov(Param::GAMMA1, Param::GAMMA2));
    set<E1_R_COV>(cov(Param::GAMMA1, Param::KAPPA));
    set<E1_N_COV>(cov(Param::GAMMA1, Param::N));

    set<E2_E2_COV>(cov(Param::GAMMA2, Param::GAMMA2));
    set<E2_R_COV>(cov(Param::GAMMA2, Param::KAPPA));
    set<E2_N_COV>(cov(Param::GAMMA2, Param::N));

    set<R_R_COV>(radiusVar);
    set<R_N_COV>(cov(Param::KAPPA, Param::N));

    set<N_N_COV>(cov(Param::N, Param::N));
}

void multifit::SmallGalaxyModelPhotometry::defineSchema(
    Schema::Ptr schema
) {
    Base::defineSchema(schema);
    schema->add(SchemaEntry("amplitude", AMPLITUDE, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e1", E1, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("e2", E2, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("r", R, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("sersic", N, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("ampAmpCov", AMP_AMP_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("ampE1Cov", AMP_E1_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("ampE2Cov", AMP_E2_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("ampRCov", AMP_R_COV, Schema::DOUBLE, 1));
    schema->add(SchemaEntry("ampNCov", AMP_N_COV, Schema::DOUBLE, 1));
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
