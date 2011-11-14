#include "lsst/meas/multifit/ShapeletModelPhotometry.h"
#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/algorithms/Measure.h"

namespace pexPolicy = lsst::pex::policy;
namespace afwDetection = lsst::afw::detection;

namespace lsst {
namespace meas {
namespace multifit {

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
                                            _nCoeff));
}

ShapeletModelPhotometry::ShapeletModelPhotometry(
    boost::int64_t const status, int nCoeff
) : lsst::afw::detection::Photometry(), _nCoeff(nCoeff)
{
    init(); 
    set<STATUS>(status);
    set<FLUX>(std::numeric_limits<double>::quiet_NaN());
    set<FLUX_ERR>(std::numeric_limits<double>::quiet_NaN());
    set<E1>(std::numeric_limits<double>::quiet_NaN());
    set<E2>(std::numeric_limits<double>::quiet_NaN());
    set<RADIUS>(std::numeric_limits<double>::quiet_NaN());
    for(int i = 0; i < _nCoeff; ++i) {
        set<COEFFICIENTS>(i, std::numeric_limits<double>::quiet_NaN());
    }
}


ShapeletModelPhotometry::ShapeletModelPhotometry(
    boost::int64_t status,
    int nCoeff,
    double flux, double fluxErr, 
    double e1, double e2, double radius,
    ndarray::Array<double const, 1,1> coefficients
) : afw::detection::Photometry(),
    _nCoeff(nCoeff)
{
    init();

    set<STATUS>(status);
    set<FLUX>(flux);
    set<FLUX_ERR>(fluxErr);
    set<RADIUS>(radius); 
    set<E1>(e1);
    set<E2>(e2);           
    for(int i = 0; i < _nCoeff; ++i) {
        set<COEFFICIENTS>(i, coefficients[i]);
    }
}


/**
 * @brief An algorithm to perform shapelet model photometry
 */
template<typename ExposureT>
class ShapeletModelAlgorithm : public meas::algorithms::Algorithm<afwDetection::Photometry, ExposureT>
{
public:
    typedef meas::algorithms::Algorithm<afwDetection::Photometry, ExposureT> AlgorithmT;
    
    /// Ctor
    ShapeletModelAlgorithm(
        SourceMeasurement::Options const& options=SourceMeasurement::readPolicy(pex::policy::Policy()),
        int nCoeff=0
        ) : AlgorithmT(), _options(options), _nCoeff(SourceMeasurement::computeCoefficientCount(options)) {}
    
    virtual std::string getName() const { return "SHAPELET_MODEL"; }
    
    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<ShapeletModelAlgorithm<ExposureT> >(_options, _nCoeff);
    }
    
    virtual void configure(pexPolicy::Policy const& policy) {
        _options = SourceMeasurement::readPolicy(policy);
        _nCoeff = SourceMeasurement::computeCoefficientCount(_options);
    }
    
    virtual PTR(afwDetection::Photometry) measureNull(void) const {
        // XXX What is the correct flag value for "I can't measure this"???
        boost::int64_t const flag = algorithms::Flags::SHAPELET_PHOTOM_BAD;
        return boost::make_shared<ShapeletModelPhotometry>(flag, _nCoeff);
    }
    
    virtual PTR(afwDetection::Photometry) measureSingle(
        afwDetection::Source const&,
        afwDetection::Source const&,
        meas::algorithms::ExposurePatch<ExposureT> const&
        ) const;
private:
    SourceMeasurement::Options _options;
    int _nCoeff;
};


template <typename ExposureT>
PTR(afw::detection::Photometry) ShapeletModelAlgorithm<ExposureT>::measureSingle( 
    afw::detection::Source const& target,
    afw::detection::Source const& source,
    meas::algorithms::ExposurePatch<ExposureT> const& patch
    ) const
{
    SourceMeasurement measurement(_options);

    int status = measurement.measure(patch.getExposure(), source); 
    if (status & algorithms::Flags::SHAPELET_PHOTOM_BAD) {
        return ShapeletModelPhotometry::Ptr(
            new ShapeletModelPhotometry(status, _nCoeff)
        );
    }

    EllipseCore const & core(measurement.getEllipse().getCore());
    return ShapeletModelPhotometry::Ptr(
        new ShapeletModelPhotometry(
            status, _nCoeff,
            measurement.getFlux(), measurement.getFluxErr(),
            core.getE1(), core.getE2(), static_cast<double>(core.getRadius()),
            measurement.getCoefficients()
        )
    );
}

LSST_DECLARE_ALGORITHM(ShapeletModelAlgorithm, afwDetection::Photometry);

// \endcond

}}}
