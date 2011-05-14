#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/GaussNewtonOptimizer.h"
#include "lsst/meas/multifit/SimpleInterpreter.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/utils/Utils.h"
#include "lsst/utils/ieee.h"

namespace lsst { namespace meas { namespace multifit {

template <int nCoeff>
lsst::afw::image::MaskPixel lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::bitmask;

template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::isEllipticityActive;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::isRadiusActive;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::isPositionActive;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::retryWithSvd;

template <int nCoeff>
double lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::ftol;
template <int nCoeff>
double lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::gtol;
template <int nCoeff>
double lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::minStep;
template <int nCoeff>
double lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::tau;

template <int nCoeff>
int lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::maxIter;
template <int nCoeff>
int lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::nGrowFp;

template <int nCoeff>
ModelBasis::Ptr lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::basis;

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

template <int nCoeff>
ShapeletModelPhotometry<nCoeff>::ShapeletModelPhotometry(
    BaseInterpreter::ConstPtr const & interpreter
) : afw::detection::Photometry() {
    init();
    if(!interpreter) {
        return;
    }
    set<FLUX>(interpreter->computeFluxMean(0, 0));
    set<FLUX_ERR>(sqrt(interpreter->computeFluxVariance(0,0)));
    afw::geom::ellipses::Ellipse ellipse = interpreter->computeEllipseMean(0);
    EllipseCore core(ellipse.getCore());
    set<RADIUS>(static_cast<double>(core.getRadius()));
    set<E1>(core.getE1());
    set<E2>(core.getE2());        
    Eigen::VectorXd coeff = interpreter->computeCoefficientMean();
    for(int i = 0; i < nCoeff; ++i){
        set<COEFFICIENTS>(i, coeff[i]);
    }
}

afw::geom::ellipses::Ellipse makeEllipse(
    afw::detection::Source const & source, 
    afw::detection::Footprint const & fp
) {
    afw::geom::Point2D center(source.getXAstrom(), source.getYAstrom());
    if(!utils::isfinite(source.getIxx()) 
       || !utils::isfinite(source.getIyy()) 
       || !utils::isfinite(source.getIxy())
    ) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Input source moments are not finite."
        );
    }
    afw::geom::ellipses::Quadrupole quad(
        source.getIxx(), source.getIyy(), source.getIxy() 
    );

    afw::geom::ellipses::Axes axes(quad);
    if(axes.getA()<= std::numeric_limits<double>::epsilon())
        axes.setA(std::numeric_limits<double>::epsilon());
    if(axes.getB()<= std::numeric_limits<double>::epsilon())
        axes.setB(std::numeric_limits<double>::epsilon());

    double maxRadius = 2*sqrt(fp.getArea());
    if(axes.getA() > maxRadius || axes.getB() > maxRadius) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Input source radius is unreasonably large."
        );
    }
    return afw::geom::ellipses::Ellipse(axes, center);
}

afw::image::MaskPixel makeBitMask(
    pex::policy::Policy::StringArray const& maskPlaneNames 
) {
    afw::image::MaskPixel bitmask = 0;
    for (unsigned i =0; i < maskPlaneNames.size(); ++i){
        bitmask |= afw::image::Mask<afw::image::MaskPixel>::getPlaneBitMask(
            maskPlaneNames[i]
        );
    }
    return bitmask;
}

std::string makeBasisPath(int nCoeff) {
    std::string file;
    if(nCoeff == 2)
        file = "ed+00:0000.boost";
    else if(nCoeff == 8)
        file = "ed+06:2000.boost";
    else if(nCoeff == 17)
        file = "ed+15:4000.boost";

    fs::path path(utils::eups::productDir("meas_multifit"));
    path /= fs::path("data/"+file);
    return path.native_file_string();
}

template <int nCoeff>
template <typename ExposureT>
afw::detection::Photometry::Ptr ShapeletModelPhotometry<nCoeff>::doMeasure(
    CONST_PTR(ExposureT) im,
    CONST_PTR(afw::detection::Peak) peak,
    CONST_PTR(afw::detection::Source) source
) {
    if (!source) {
        return ShapeletModelPhotometry::Ptr(new ShapeletModelPhotometry()); 
    }

    afw::detection::Footprint::Ptr fp = afw::detection::growFootprint(
        *source->getFootprint(), nGrowFp
    );

    afw::geom::ellipses::Ellipse ellipse = makeEllipse(*source, *fp);

    Definition definition = Definition::make(
        *im, fp, basis, ellipse, 
        isEllipticityActive,
        isRadiusActive,
        isPositionActive,
        bitmask
    );
    Evaluator::Ptr evaluator = Evaluator::make(definition);
    GaussNewtonOptimizer optimizer;
    SimpleDistribution::Ptr distribution = optimizer.solve(
        evaluator,
        ftol, gtol, minStep, maxIter, 
        tau, retryWithSvd);

    UnifiedSimpleInterpreter::Ptr interpreter = UnifiedSimpleInterpreter::make(
        distribution, evaluator->getGrid()
    );
    return boost::make_shared<ShapeletModelPhotometry<nCoeff> >(interpreter);
}


template <int nCoeff>
bool ShapeletModelPhotometry<nCoeff>::doConfigure(
    lsst::pex::policy::Policy const & policy
) {   
    lsst::pex::policy::Policy local(policy);
    lsst::pex::policy::Policy dict;
    lsst::pex::policy::DefaultPolicyFile file("meas_multifit", "ShapeletModelPhotometryDict.paf", "policy");
    file.load(dict);

    local.mergeDefaults(dict);

    nGrowFp = local.getInt("nGrowFp");
    isEllipticityActive = local.getBool("isEllipticityActive");
    isRadiusActive = local.getBool("isRadiusActive");
    isPositionActive = local.getBool("isPositionActive");
    bitmask = makeBitMask(local.getStringArray("maskPlaneName"));

    ftol = local.getDouble("ftol");
    gtol = local.getDouble("ftol");
    minStep = local.getDouble("minStep");
    maxIter = local.getInt("maxIter");
    tau = local.getDouble("tau");
    retryWithSvd = local.getBool("retryWithSvd");

    basis = CompoundShapeletModelBasis::load(makeBasisPath(nCoeff));
    ShapeletModelBasis::setPsfShapeletOrder(local.getInt("psfShapeletOrder"));
    return true;
}

/*
 * Declare the existence of a "SHAPELET_MODEL_N" algorithm to MeasurePhotometry
 *
 * 
 *
 * @cond
 */
#define INSTANTIATE(NAME, TYPE, N_COEFF) \
    lsst::meas::algorithms::MeasurePhotometry<lsst::afw::image::Exposure<TYPE> >::declare(NAME, \
        &ShapeletModelPhotometry<N_COEFF>::doMeasure<lsst::afw::image::Exposure<TYPE> >, \
        &ShapeletModelPhotometry<N_COEFF>::doConfigure               \
    )

volatile bool isInstance[] = {
    INSTANTIATE("SHAPELET_MODEL_2", float, 2),
    INSTANTIATE("SHAPELET_MODEL_2", double, 2),
    INSTANTIATE("SHAPELET_MODEL_8", float, 8),
    INSTANTIATE("SHAPELET_MODEL_8", double, 8),
    INSTANTIATE("SHAPELET_MODEL_17", float, 17),
    INSTANTIATE("SHAPELET_MODEL_17", double, 17)
};

// \endcond

}}}
