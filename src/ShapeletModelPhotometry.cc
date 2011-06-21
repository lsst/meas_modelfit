#include "lsst/meas/multifit/ShapeletModelPhotometry.h"
#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/algorithms/Measure.h"

namespace lsst {
namespace meas {
namespace multifit {



lsst::afw::image::MaskPixel ShapeletModelPhotometry::bitmask;
ModelBasis::Ptr ShapeletModelPhotometry::basis;

bool ShapeletModelPhotometry::usePixelWeights;
bool ShapeletModelPhotometry::fitDeltaFunction;
bool ShapeletModelPhotometry::isEllipticityActive;
bool ShapeletModelPhotometry::isRadiusActive;
bool ShapeletModelPhotometry::isPositionActive;

int ShapeletModelPhotometry::psfShapeletOrder;
int ShapeletModelPhotometry::nCoeff;
int ShapeletModelPhotometry::nGrowFp;
int ShapeletModelPhotometry::nTestPoints;


void ShapeletModelPhotometry::defineSchema(lsst::afw::detection::Schema::Ptr schema) {
    schema->clear();
    schema->add(afw::detection::SchemaEntry("status",  STATUS,   afw::detection::Schema::LONG, 1));
    schema->add(afw::detection::SchemaEntry("flux",    FLUX,     afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("fluxErr", FLUX_ERR, afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("e1",      E1,       afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("e2",      E2,       afw::detection::Schema::DOUBLE, 1));
    schema->add(afw::detection::SchemaEntry("radius",  RADIUS,   afw::detection::Schema::DOUBLE, 1,
                                            "pixels"));
    schema->add(afw::detection::SchemaEntry("coefficients", COEFFICIENTS, afw::detection::Schema::DOUBLE,
                                            nCoeff));

}

ShapeletModelPhotometry::ShapeletModelPhotometry(
    boost::int64_t const status
) : lsst::afw::detection::Photometry(),
    _flag(status)
{
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


ShapeletModelPhotometry::ShapeletModelPhotometry(
    boost::int64_t status,  
    double flux, double fluxErr, 
    double e1, double e2, double radius,
    ndarray::Array<double const, 1,1> coeff
) : afw::detection::Photometry(),
    _flag(status)
{
    init();

    set<STATUS>(status);
    set<FLUX>(flux);
    set<FLUX_ERR>(fluxErr);
    set<RADIUS>(radius); 
    set<E1>(e1);
    set<E2>(e2);           
    for(int i = 0; i < nCoeff; ++i) {
        set<COEFFICIENTS>(i, coeff[i]);
    }
}

template <typename ExposureT>
ShapeletModelPhotometry::Photometry::Ptr ShapeletModelPhotometry::doMeasure(
    CONST_PTR(ExposureT) exp,
    CONST_PTR(afw::detection::Peak) peak,
    CONST_PTR(afw::detection::Source) source
) {
    SourceMeasurement measurement(basis, psfShapeletOrder, 
        nTestPoints, nGrowFp, usePixelWeights, fitDeltaFunction,
        isEllipticityActive, isRadiusActive, isPositionActive,
        bitmask
    );

    int status = measurement.measure(exp, source); 
    if (status & algorithms::Flags::SHAPELET_PHOTOM_BAD) {
        return ShapeletModelPhotometry::Ptr(
            new ShapeletModelPhotometry(status)
        );
    }

    afw::geom::ellipses::Ellipse::ConstPtr ellipse(measurement.getEllipse());
    EllipseCore const & core(ellipse->getCore());
    return ShapeletModelPhotometry::Ptr(
        new ShapeletModelPhotometry(            
            status,
            measurement.getFlux(), measurement.getFluxErr(),
            core.getE1(), core.getE2(), static_cast<double>(core.getRadius()),
            measurement.getCoeff()
        )
    );
}

bool ShapeletModelPhotometry::doConfigure(
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
    bitmask = lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask(
        local.getStringArray("maskPlaneName")
    );
    fitDeltaFunction = local.getBool("fitDeltaFunction");

    int basisSize = local.getInt("basisSize");
    basis = SourceMeasurement::loadBasis(basisSize);
    nCoeff = fitDeltaFunction ? basisSize + 1: basisSize;

    usePixelWeights = local.getBool("usePixelWeights");
    
    nTestPoints = local.getInt("nTestPoints");
    psfShapeletOrder = local.getInt("psfShapeletOrder");
    
    return true;
}



/*
 * Declare the existence of a "SHAPELET_MODEL" algorithm to MeasurePhotometry
 *
 * @cond
 */
#define INSTANTIATE(NAME, TYPE) \
    lsst::meas::algorithms::MeasurePhotometry<lsst::afw::image::Exposure<TYPE> >::declare(NAME, \
        &ShapeletModelPhotometry::doMeasure<lsst::afw::image::Exposure<TYPE> >, \
        &ShapeletModelPhotometry::doConfigure               \
    )

volatile bool isInstance[] = {
    INSTANTIATE("SHAPELET_MODEL", float),
    INSTANTIATE("SHAPELET_MODEL", double),
};

// \endcond

}}}
