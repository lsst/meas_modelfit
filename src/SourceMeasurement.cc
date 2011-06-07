#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/meas/multifit/SimpleInterpreter.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/utils/Utils.h"
#include "lsst/utils/ieee.h"

namespace lsst { namespace meas { namespace multifit {


template <int basisSize>
bool lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::usePixelWeights;
template <int basisSize>
lsst::afw::image::MaskPixel lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::bitmask;
template <int basisSize>
ModelBasis::Ptr lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::basis;
template <int basisSize>
int lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::nCoeff;
template <int basisSize>
bool lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::fitDeltaFunction;
template <int basisSize>
int lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::nGrowFp;
template <int basisSize>
bool lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::isEllipticityActive;
template <int basisSize>
bool lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::isRadiusActive;
template <int basisSize>
bool lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::isPositionActive;


template <int basisSize>
int lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::nTestPoints;

#if 0 
//these are needed for the GaussNewtonoptimizer
template <int basisSize>
bool lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::retryWithSvd;
template <int basisSize>
double lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::ftol;
template <int basisSize>
double lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::gtol;
template <int basisSize>
double lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::minStep;
template <int basisSize>
double lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::tau;
template <int basisSize>
int lsst::meas::multifit::ShapeletModelPhotometry<basisSize>::maxIter;
#endif




template <int basisSize>
void ShapeletModelPhotometry<basisSize>::defineSchema(lsst::afw::detection::Schema::Ptr schema) {
    schema->clear();
    schema->add(afw::detection::SchemaEntry("status",  STATUS,   afw::detection::Schema::INT, 1));
    schema->add(afw::detection::SchemaEntry("flux",    FLUX,     afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("fluxErr", FLUX_ERR, afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("e1",      E1,       afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("e2",      E2,       afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("radius",  RADIUS,   afw::detection::Schema::DOUBLE, 1,
                                            "pixels"));
    schema->add(afw::detection::SchemaEntry("coefficients", COEFFICIENTS, afw::detection::Schema::DOUBLE,
                                            nCoeff));

}

template <int basisSize>
ShapeletModelPhotometry<basisSize>::ShapeletModelPhotometry(
    int const status
) : lsst::afw::detection::Photometry() {
    init(); 
    set<STATUS>(status);
    set<FLUX>(std::numeric_limits<double>::quiet_NaN());
    set<FLUX_ERR>(std::numeric_limits<double>::quiet_NaN());
    set<E1>(std::numeric_limits<double>::quiet_NaN());
    set<E2>(std::numeric_limits<double>::quiet_NaN());
    set<RADIUS>(std::numeric_limits<double>::quiet_NaN());
    for(int i = 0; i < nCoeff; ++i) {
        set<COEFFICIENTS>(i, std::numeric_limits<double>::quiet_NaN());
    }
}

template <int basisSize>
ShapeletModelPhotometry<basisSize>::ShapeletModelPhotometry(
    Evaluator::Ptr const & evaluator,
    ndarray::Array<const double, 1, 1> const & param,
    ndarray::Array<const double, 1, 1> const & coeff,
    ndarray::Array<const double, 2, 1> const & covar,
    int const status
) : afw::detection::Photometry() {
    init();
    
    set<STATUS>(status);
    if(status & OPTIMIZER_FAILED) {
        set<FLUX>(std::numeric_limits<double>::quiet_NaN());
        set<FLUX_ERR>(std::numeric_limits<double>::quiet_NaN());
        set<E1>(std::numeric_limits<double>::quiet_NaN());
        set<E2>(std::numeric_limits<double>::quiet_NaN());
        set<RADIUS>(std::numeric_limits<double>::quiet_NaN());
        for(int i = 0; i < nCoeff; ++i) {
            set<COEFFICIENTS>(i, std::numeric_limits<double>::quiet_NaN());
        }        
        return;
    }
    int starId=0; 
    double flux=0, fluxVar=0;

    if((status & GALAXY_MODEL_FAILED) == 0) {
        starId = 1;
        grid::Object const & object = evaluator->getGrid()->objects[0];
        grid::Source const & source = evaluator->getGrid()->sources[0];
        flux += source.computeFluxMean(param, coeff);
        fluxVar += source.computeFluxVariance(param, covar);
        afw::geom::ellipses::Ellipse ellipse = object.makeEllipse(param.begin());
        EllipseCore core(ellipse.getCore());
        set<RADIUS>(static_cast<double>(core.getRadius())); 
        set<E1>(core.getE1());
        set<E2>(core.getE2());       

        for(int i = 0; i < basisSize; ++i) {
            set<COEFFICIENTS>(i, coeff[i]);
        }
    } else {
        set<RADIUS>(0.); 
        set<E1>(0.);
        set<E2>(0.);   
        for(int i = 0; i < basisSize; ++i) {
            set<COEFFICIENTS>(i, std::numeric_limits<double>::quiet_NaN());
        }
    }
    if(fitDeltaFunction) {  
        grid::Source const & star = evaluator->getGrid()->sources[starId];
        flux += star.computeFluxMean(param, coeff);
        fluxVar += star.computeFluxVariance(param, covar);        
        set<COEFFICIENTS>(basisSize, coeff[star.getCoefficientOffset()]);
    }
    set<FLUX>(flux);
    set<FLUX_ERR>(sqrt(fluxVar));

}

#if 0 
template <int basisSize>
ShapeletModelPhotometry<basisSize>::ShapeletModelPhotometry(
    GaussNewtonOptimizer & optimizer,
    BaseEvaluator::Ptr const &
) : afw::detection::Photometry() {
    init();
    SimpleDistribution::Ptr distribution = optimizer.solve(
        evaluator,
        ftol, gtol, minStep, maxIter, 
        tau, retryWithSvd);

    UnifiedSimpleInterpreter::Ptr interpreter = UnifiedSimpleInterpreter::make(
        distribution, evaluator->getGrid()
    );
    int status = (optimizer.didConverge())? 0 : OPTIMIZER_FAILED;

    set<STATUS>(status);
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
#endif

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

std::string makeBasisPath(int basisSize) {
    std::string file;
    if(basisSize == 2)
        file = "ed+00:0000.boost";
    else if(basisSize == 8)
        file = "ed+06:2000.boost";
    else if(basisSize == 17)
        file = "ed+15:4000.boost";
    else {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Unsupported number of coefficients. No corresponding basis exists on file"
        );
    }
    fs::path path(utils::eups::productDir("meas_multifit"));
    path /= fs::path("data/"+file);
    return path.native_file_string();
}

template <int basisSize>
template <typename ExposureT>
afw::detection::Photometry::Ptr ShapeletModelPhotometry<basisSize>::doMeasure(
    CONST_PTR(ExposureT) im,
    CONST_PTR(afw::detection::Peak) peak,
    CONST_PTR(afw::detection::Source) source
) {
  std::cerr << "starting source " << source->getSourceId() << std::endl;
    if (!source) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(NO_SOURCE));
    }
    if (!source->getFootprint()) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(NO_FOOTPRINT));
    }
    if (!basis) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(NO_BASIS));
    }
    if (!im) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(NO_EXPOSURE));
    }
    if (!im->getPsf()) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(NO_PSF));
    }

#if 0
    if (source->getSourceId() != 2071) {
        // Just to see if it's specifically this object that matters.
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(NO_EXPOSURE));
    }
#endif

    afw::detection::Footprint::Ptr fp = afw::detection::growFootprint(
        *source->getFootprint(), nGrowFp
    );

    boost::scoped_ptr<afw::geom::ellipses::Ellipse> ellipse;

    try{
        ellipse.reset(new Ellipse(makeEllipse(*source, *fp)));
    } catch(lsst::pex::exceptions::InvalidParameterException e) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(BAD_INITIAL_MOMENTS));
    }

    
    Definition definition;
    definition.frames.insert(definition::Frame::make(0, *im, fp, bitmask, usePixelWeights));
    //fit both a point source and galaxy model for the same object
    definition::Object galaxy = definition::Object::makeGalaxy(
        0, basis, *ellipse, 
        isEllipticityActive, 
        isRadiusActive, 
        isPositionActive
    );
    definition::Object star(1);
    //link the position object of the two models
    star.getPosition() = galaxy.getPosition();

    definition.objects.insert(galaxy);
    definition.objects.insert(star);


    ndarray::Array<const double, 1, 1> param, coeff;
    ndarray::Array<const double, 2, 1> covar;
    int status = 0;
    Evaluator::Ptr evaluator = Evaluator::make(definition);
    BruteForceSourceOptimizer optimizer;
    bool success = optimizer.solve(evaluator, nTestPoints);
    if(success) {
        param = optimizer.getBestParameters();
        coeff = optimizer.getBestCoefficients();
        covar = optimizer.getCoefficientCovariance();
    }
    else if (fitDeltaFunction) { 
        status |= GALAXY_MODEL_FAILED;

        //try fitting just the point source model
        definition.objects.erase(galaxy.id);
        evaluator = Evaluator::make(definition);
        Evaluation evaluation(evaluator);
        try{
            param = evaluation.getParameters();
            coeff = evaluation.getCoefficients();
            covar = evaluation.getCoefficientFisherMatrix();
        }
        catch (...){
            status |= static_cast<int>(OPTIMIZER_FAILED);
        }
    }
    return ShapeletModelPhotometry<basisSize>::Ptr(
        new ShapeletModelPhotometry(evaluator, param, coeff, covar, status)
    ); 
}


template <int basisSize>
bool ShapeletModelPhotometry<basisSize>::doConfigure(
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
    fitDeltaFunction = local.getBool("fitDeltaFunction");
    nCoeff = fitDeltaFunction ? basisSize + 1: basisSize;

    usePixelWeights = local.getBool("usePixelWeights");
    
    nTestPoints = local.getInt("nTestPoints");
    
#if 0
    ftol = local.getDouble("ftol");
    gtol = local.getDouble("ftol");
    minStep = local.getDouble("minStep");
    maxIter = local.getInt("maxIter");
    tau = local.getDouble("tau");
    retryWithSvd = local.getBool("retryWithSvd");
#endif

    basis = CompoundShapeletModelBasis::load(makeBasisPath(basisSize));
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

template class lsst::meas::multifit::ShapeletModelPhotometry<2>;
template class lsst::meas::multifit::ShapeletModelPhotometry<8>;
template class lsst::meas::multifit::ShapeletModelPhotometry<17>;


// \endcond

}}}
