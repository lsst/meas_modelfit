#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/utils/Utils.h"
#include "lsst/utils/ieee.h"

#include "Eigen/Cholesky"

#define LSST_MAX_DEBUG 0
#include "lsst/pex/logging/Debug.h"
#include <iostream>

namespace lsst { namespace meas { namespace multifit {

namespace {

static double const SQRT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());

afw::geom::ellipses::Quadrupole naiveDeconvolve(
    afw::geom::ellipses::BaseCore const & ellipse,
    afw::geom::ellipses::BaseCore const & psf
) {
    afw::geom::ellipses::Quadrupole moments(ellipse);
    afw::geom::ellipses::Quadrupole psfMoments(psf);
    double ixx = moments.getIXX() - psfMoments.getIXX();
    double iyy = moments.getIYY() - psfMoments.getIYY();
    double ixy = moments.getIXY() - psfMoments.getIXY();
    if ((ixx < SQRT_EPS) || (iyy < SQRT_EPS) || (ixy * ixy > ixx * iyy)) {
        ixx = SQRT_EPS;
        iyy = SQRT_EPS;
        ixy = 0.0;
    }
    return afw::geom::ellipses::Quadrupole(ixx, iyy, ixy);
}

} // anonymous

SourceMeasurement::Options SourceMeasurement::readPolicy(
    lsst::pex::policy::Policy const & policy
) {   
    lsst::pex::policy::Policy local(policy);
    lsst::pex::policy::Policy dict;
    lsst::pex::policy::DefaultPolicyFile file("meas_multifit", "ShapeletModelPhotometryDict.paf", "policy");
    file.load(dict);

    local.mergeDefaults(dict);

    Options options;
    options.fitDeltaFunction = local.getBool("fitDeltaFunction");
    options.fitExponential = local.getBool("fitExponential");
    options.fitDeVaucouleur = local.getBool("fitDeVaucouleur");
    options.shapeletOrder = local.getInt("shapeletOrder");
    options.psfShapeletOrder = local.getInt("psfShapeletOrder");
    options.nGrowFp = local.getInt("nGrowFp");
    options.usePixelWeights = local.getBool("usePixelWeights");
    options.ellipticityStepSize = local.getDouble("ellipticityStepSize");
    options.ellipticityStepCount = local.getInt("ellipticityStepCount");
    options.radiusMaxFactor = local.getDouble("radiusMaxFactor");
    options.radiusMinFactor = local.getDouble("radiusMinFactor");
    options.radiusStepCount = local.getInt("radiusStepCount");
    options.maskPlaneNames = local.getStringArray("maskPlaneName");

    if (options.fitExponential && !SourceMeasurement::getExponentialBasis()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Could not load exponential model from meas_multifit/data.");
    }
    if (options.fitDeVaucouleur && !SourceMeasurement::getDeVaucouleurBasis()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Could not load de Vaucouleur model from meas_multifit/data.");
    }
    return options;
}

afw::geom::ellipses::Ellipse SourceMeasurement::makeEllipse(
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

SourceMeasurement::SourceMeasurement(Options const & options) : 
    _options(options),
    _bitmask(afw::image::Mask<afw::image::MaskPixel>::getPlaneBitMask(options.maskPlaneNames)),
    _ellipse(EllipseCore())
{
    int coefficientSize = computeCoefficientSize(options);
    _integration = ndarray::allocate(coefficientSize);
    _coefficients = ndarray::allocate(coefficientSize);
    _covariance = ndarray::allocate(coefficientSize, coefficientSize);
    _parameters = ndarray::allocate(3);
    _integration.deep() = 0.0;
    int offset = 0;
    if (_options.fitDeltaFunction) {
        _integration[offset] = 1.0;
        ++offset;
    }
    if (_options.fitExponential) {
        getExponentialBasis()->integrate(_integration[ndarray::view(offset, offset+1)]);
        ++offset;
    }
    if (_options.fitDeVaucouleur) {
        getDeVaucouleurBasis()->integrate(_integration[ndarray::view(offset, offset+1)]);
        ++offset;
    }
    if (_options.shapeletOrder >= 0) {
        CompoundShapeletBuilder::ComponentVector components;
        components.push_back(ShapeletModelBasis::make(_options.shapeletOrder));
        CompoundShapeletBuilder builder(components);
        builder.normalizeFlux(0);
        ndarray::Array<Pixel,2,2> constraintMatrix(ndarray::allocate(2, builder.getSize()));
        constraintMatrix.deep() = 0.0;
        builder.integrate(constraintMatrix[1]);
        std::swap(constraintMatrix[1][0], constraintMatrix[0][0]);
        ndarray::Array<Pixel,1,1> constraintVector(ndarray::allocate(2));
        constraintVector.deep() = 0.0;
        builder.setConstraint(constraintMatrix, constraintVector);
        _shapeletBasis = builder.build();
        _shapeletBasis->integrate(_integration[ndarray::view(offset, offset + _shapeletBasis->getSize())]);
    }
}

CompoundShapeletModelBasis::Ptr SourceMeasurement::loadBasis(std::string const & name) {
    fs::path path(utils::eups::productDir("meas_multifit"));
    path /= fs::path("data");
    path /= name + ".boost";
    return CompoundShapeletModelBasis::load(path.native_file_string());    
}

CompoundShapeletModelBasis::Ptr SourceMeasurement::getExponentialBasis() {
    static CompoundShapeletModelBasis::Ptr cached = loadBasis("exponential");
    return cached;
}

CompoundShapeletModelBasis::Ptr SourceMeasurement::getDeVaucouleurBasis() {
    static CompoundShapeletModelBasis::Ptr cached = loadBasis("deVaucouleur");
    return cached;
}

void SourceMeasurement::addObjectsToDefinition(
    Definition & def, lsst::afw::geom::ellipses::Ellipse const & ellipse
) const {
    definition::PositionComponent::Ptr position = definition::PositionComponent::make(
        ellipse.getCenter(), false
    );
    EllipseCore ellipseCore(ellipse.getCore());
    definition::RadiusComponent::Ptr radius = definition::RadiusComponent::make(
        ellipseCore.getRadius(), true
    );
    definition::EllipticityComponent::Ptr ellipticity = definition::EllipticityComponent::make(
        ellipseCore.getEllipticity(), true
    );
    if (_options.fitDeltaFunction) {
        definition::ObjectComponent obj(DELTAFUNCTION_ID);
        obj.getPosition() = position;
        def.objects.insert(obj);
    }
    if (_options.fitExponential) {
        definition::ObjectComponent obj(EXPONENTIAL_ID);
        obj.getPosition() = position;
        obj.getRadius() = radius;
        obj.getEllipticity() = ellipticity;
        obj.setBasis(getExponentialBasis());
        def.objects.insert(obj);

    }
    if (_options.fitDeVaucouleur) {
        definition::ObjectComponent obj(DEVAUCOULEUR_ID);
        obj.getPosition() = position;
        obj.getRadius() = radius;
        obj.getEllipticity() = ellipticity;
        obj.setBasis(getDeVaucouleurBasis());
        def.objects.insert(obj);
    }
    if (_options.shapeletOrder >= 0) {
        definition::ObjectComponent obj(SHAPELET_ID);
        obj.getPosition() = position;
        obj.getRadius() = radius;
        obj.getEllipticity() = ellipticity;
        obj.setBasis(_shapeletBasis);
        def.objects.insert(obj);
    }
}

void SourceMeasurement::solve(double e1, double e2, double radius, double & best) {
    pex::logging::Debug log("photometry.multifit", LSST_MAX_DEBUG);
    assert(_evaluator->getParameterSize() == 3);
    log.debug(4, boost::format("Testing %f %f %f.") % e1 % e2 % radius);
    ndarray::Array<double,1,1> parameters(ndarray::allocate(_evaluator->getParameterSize()));
    parameters[_evaluator->getGrid()->radii[0].offset] = radius;
    parameters[_evaluator->getGrid()->ellipticities[0].offset] = e1;
    parameters[_evaluator->getGrid()->ellipticities[0].offset + 1] = e2;
    Evaluation evaluation(_evaluator, parameters);
    double objective;
    try {
        objective = evaluation.getObjectiveValue();
        log.debug(6, boost::format("Objective value is %f") % objective);
    } catch (...) {
        return;
    }
    if (objective < best) {
        _status &= ~algorithms::Flags::SHAPELET_PHOTOM_GALAXY_FAIL;
        //double flux = grid::SourceComponent::computeFlux(_integration, evaluation.getCoefficients());
        //double condition = flux / ndarray::viewAsEigen(evaluation.getCoefficients()).norm();
        //if (condition < 1E-10) {
        //    if (!(_status & algorithms::Flags::SHAPELET_PHOTOM_INVERSION_UNSAFE)) return;
        //} else {
        //    _status &= ~algorithms::Flags::SHAPELET_PHOTOM_INVERSION_UNSAFE;
        // }
        _parameters.deep() = parameters;
        _coefficients.deep() = evaluation.getCoefficients();
        _covariance.deep() = evaluation.getCoefficientFisherMatrix();
        best = objective;
        _ellipse.getCore() = EllipseCore(e1, e2, radius);
    }
}

void SourceMeasurement::optimize(Ellipse const & initialEllipse) {
    EllipseCore maxEllipse(initialEllipse.getCore());
    EllipseCore minEllipse(maxEllipse);
    afw::geom::Ellipse psfEllipse = _evaluator->getGrid()->sources[0].getLocalPsf()->computeMoments();
    maxEllipse.scale(_options.radiusMaxFactor);
    minEllipse.scale(_options.radiusMinFactor);
    maxEllipse = naiveDeconvolve(maxEllipse, psfEllipse.getCore());
    minEllipse = naiveDeconvolve(minEllipse, psfEllipse.getCore());
    EllipseCore::ParameterVector base(minEllipse.getParameterVector());
    EllipseCore::ParameterVector delta = maxEllipse.getParameterVector() - minEllipse.getParameterVector();
    delta /= _options.radiusStepCount;
    double best = std::numeric_limits<double>::infinity();
    _status |= algorithms::Flags::SHAPELET_PHOTOM_GALAXY_FAIL;
    for (int iR = 0; iR < _options.radiusStepCount; ++iR, base += delta) {
        if (base[EllipseCore::RADIUS] < SQRT_EPS) {
            solve(0.0, 0.0, SQRT_EPS, best);
            continue;
        }
        for (int iE1 = -_options.ellipticityStepCount; iE1 <= _options.ellipticityStepCount; ++iE1) {
            for (int iE2 = -_options.ellipticityStepCount; iE2 <= _options.ellipticityStepCount; ++iE2) {
                solve(
                    base[EllipseCore::E1] + iE1 * _options.ellipticityStepSize,
                    base[EllipseCore::E2] + iE2 * _options.ellipticityStepSize,
                    base[EllipseCore::RADIUS],
                    best
                );
            }
        }
    }
    if (_status & algorithms::Flags::SHAPELET_PHOTOM_GALAXY_FAIL) {
        _ellipse = initialEllipse;
        _flux = std::numeric_limits<double>::quiet_NaN();
        _fluxErr = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    Eigen::LDLT<Eigen::MatrixXd> ldlt(ndarray::viewAsTransposedEigen(_covariance));
    ndarray::viewAsTransposedEigen(_covariance).setIdentity();
    ldlt.solveInPlace(ndarray::viewAsTransposedEigen(_covariance).setIdentity());
    Grid::Ptr grid = _evaluator->getGrid();
    _flux = grid::SourceComponent::computeFlux(_integration, _coefficients);
    _fluxErr = std::sqrt(grid::SourceComponent::computeFluxVariance(_integration, _covariance));
}

template <typename ExposureT>
int SourceMeasurement::measure(
    CONST_PTR(ExposureT) exp,
    CONST_PTR(afw::detection::Source) source
) {
    pex::logging::Debug log("photometry.multifit", LSST_MAX_DEBUG);
    _status = 0;    
    if (!source) {
        _status |= algorithms::Flags::PHOTOM_NO_SOURCE;
        return _status;
    }
    log.debug(1, boost::format("Processing source %lld") % source->getSourceId());
    if (!source->getFootprint()) {
        _status |= algorithms::Flags::PHOTOM_NO_FOOTPRINT;
        return _status;
    }
    if (!exp->getPsf()) {
        _status |= algorithms::Flags::PHOTOM_NO_PSF;
        return _status;
    }
    ShapeletModelBasis::setPsfShapeletOrder(_options.psfShapeletOrder);
    afw::detection::Footprint::Ptr fp = afw::detection::growFootprint(
        *source->getFootprint(), _options.nGrowFp
    );
    boost::scoped_ptr<afw::geom::ellipses::Ellipse> ellipse;
    try{
        ellipse.reset(new Ellipse(makeEllipse(*source, *fp)));
        _ellipse.setCenter(ellipse->getCenter());
    } catch(lsst::pex::exceptions::InvalidParameterException e) {
        _status |= algorithms::Flags::SHAPELET_PHOTOM_BAD_MOMENTS;
        return _status;
    }
    Definition def;
    def.frames.insert(definition::Frame::make(0, *exp, fp, _bitmask, _options.usePixelWeights));
    _fp = def.frames[0].getFootprint();
    if (_fp->getArea() == 0) {
        _status |= algorithms::Flags::PHOTOM_NO_FOOTPRINT;
        return _status;
    }
    addObjectsToDefinition(def, *ellipse);
    _evaluator = Evaluator::make(def);

    optimize(*ellipse);

    return _status;

}
template int SourceMeasurement::measure(
    CONST_PTR(afw::image::Exposure<float>) exp,
    CONST_PTR(afw::detection::Source) source
);

template int SourceMeasurement::measure(
    CONST_PTR(afw::image::Exposure<double>) exp,
    CONST_PTR(afw::detection::Source) source
);

}}}
