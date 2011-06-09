#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/utils/Utils.h"
#include "lsst/utils/ieee.h"

namespace lsst { namespace meas { namespace multifit {


template <int nCoeff>
lsst::afw::image::MaskPixel lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::bitmask;
template <int nCoeff>
ModelBasis::Ptr lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::basis;
template <int nCoeff>
int lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::nGrowFp;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::isEllipticityActive;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::isRadiusActive;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::isPositionActive;
template <int nCoeff>
bool lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::addPointSource;


template <int nCoeff>
int lsst::meas::multifit::ShapeletModelPhotometry<nCoeff>::nTestPoints;

template <int nCoeff>
void ShapeletModelPhotometry<nCoeff>::defineSchema(lsst::afw::detection::Schema::Ptr schema) {
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

template <int nCoeff>
ShapeletModelPhotometry<nCoeff>::ShapeletModelPhotometry(
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

template <int nCoeff>
template <typename ExposureT>
afw::detection::Photometry::Ptr ShapeletModelPhotometry<nCoeff>::doMeasure(
    CONST_PTR(ExposureT) im,
    CONST_PTR(afw::detection::Peak) peak,
    CONST_PTR(afw::detection::Source) source
) {
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

    afw::detection::Footprint::Ptr fp = afw::detection::growFootprint(
        *source->getFootprint(), nGrowFp
    );

    boost::scoped_ptr<afw::geom::ellipses::Ellipse> ellipse;

    try{
        ellipse.reset(new Ellipse(makeEllipse(*source, *fp)));
    } catch(lsst::pex::exceptions::InvalidParameterException e) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(BAD_INITIAL_MOMENTS));
    }


    Definition def;
    def.frames.insert(definition::Frame::make(0, *im, fp, bitmask));
    def.objects.insert(
        definition::ObjectComponent::makeGalaxy(
            0, basis, *ellipse,
            isEllipticityActive,
            isRadiusActive,
            isPositionActive
        )
    );
    if (addPointSource) {
        def.objects.insert(definition::ObjectComponent(1));
        def.objects[1].getPositionElement() = def.objects[0].getPositionElement();
    }
    Evaluator::Ptr evaluator = Evaluator::make(def);
#if 0 // TODO
    BruteForceSourceOptimizer optimizer;
    bool success = optimizer.solve(evaluator, nTestPoints);
    if(!success) {
        return boost::make_shared<ShapeletModelPhotometry>(static_cast<int>(OPTIMIZER_FAILED));
    }
    return ShapeletModelPhotometry::Ptr(
        new ShapeletModelPhotometry(optimizer, evaluator)
    ); 
#endif
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
    addPointSource = local.getBool("addPointSource");
    bitmask = makeBitMask(local.getStringArray("maskPlaneName"));
    
    
    nTestPoints = local.getInt("nTestPoints");

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

template class lsst::meas::multifit::ShapeletModelPhotometry<2>;
template class lsst::meas::multifit::ShapeletModelPhotometry<8>;
template class lsst::meas::multifit::ShapeletModelPhotometry<17>;


// \endcond

}}}
