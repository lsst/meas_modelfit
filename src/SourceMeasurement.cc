#include "lsst/meas/multifit/SourceMeasurement.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/utils/Utils.h"
#include "lsst/utils/ieee.h"

#define LSST_MAX_DEBUG 0
#include "lsst/pex/logging/Debug.h"

namespace lsst { namespace meas { namespace multifit {

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

ModelBasis::Ptr SourceMeasurement::loadBasis(std::string const & name) {
    fs::path path(utils::eups::productDir("meas_multifit"));
    path /= fs::path("data");
    path /= name;
    if(!fs::exists(path)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "No corresponding basis exists on file"
        );
    }

    return CompoundShapeletModelBasis::load(path.native_file_string());
}

ModelBasis::Ptr SourceMeasurement::loadBasis(int basisSize) {
    std::string name;
    if(basisSize == 2)
        name = "ed+00:0000.boost";
    else if(basisSize == 8)
        name = "ed+06:2000.boost";
    else if(basisSize == 17)
        name = "ed+15:4000.boost";
    else {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "No corresponding basis exists on file"
        );
    }

    return loadBasis(name);
}

SourceMeasurement::SourceMeasurement(
    int basisSize, int psfShapeletOrder,
    int nTestPoints, int nGrowFp, 
    bool usePixelWeights, bool fitDeltaFunction,
    bool isEllipticityActive, 
    bool isRadiusActive, 
    bool isPositionActive,
    std::vector<std::string> const & maskPlaneNames
) : _usePixelWeights(usePixelWeights), _fitDeltaFunction(fitDeltaFunction),
    _isEllipticityActive(isEllipticityActive),
    _isRadiusActive(isRadiusActive),
    _isPositionActive(isPositionActive),
    _bitmask(afw::image::Mask<afw::image::MaskPixel>::getPlaneBitMask(maskPlaneNames)),
    _nTestPoints(nTestPoints), _nGrowFp(nGrowFp), 
    _psfShapeletOrder(psfShapeletOrder), 
    _basis(loadBasis(basisSize))
{}

SourceMeasurement::SourceMeasurement(
    ModelBasis::Ptr basis, int psfShapeletOrder,
    int nTestPoints, int nGrowFp, 
    bool usePixelWeights, bool fitDeltaFunction,
    bool isEllipticityActive, 
    bool isRadiusActive, 
    bool isPositionActive,
    lsst::afw::image::MaskPixel bitmask
) : _usePixelWeights(usePixelWeights), _fitDeltaFunction(fitDeltaFunction),
    _isEllipticityActive(isEllipticityActive),
    _isRadiusActive(isRadiusActive),
    _isPositionActive(isPositionActive),
    _bitmask(bitmask),
    _nTestPoints(nTestPoints), _nGrowFp(nGrowFp), 
    _psfShapeletOrder(psfShapeletOrder),
    _basis(basis)
{}

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
    if (!_basis) {
        _status |= algorithms::Flags::SHAPELET_PHOTOM_NO_BASIS;
        return _status;
    }
    if (!exp->getPsf()) {
        _status |= algorithms::Flags::PHOTOM_NO_PSF;
        return _status;
    }

    ShapeletModelBasis::setPsfShapeletOrder(_psfShapeletOrder);

    afw::detection::Footprint::Ptr fp = afw::detection::growFootprint(
        *source->getFootprint(), _nGrowFp
    );

    boost::scoped_ptr<afw::geom::ellipses::Ellipse> ellipse;

    try{
        ellipse.reset(new Ellipse(makeEllipse(*source, *fp)));
    } catch(lsst::pex::exceptions::InvalidParameterException e) {
        _status |= algorithms::Flags::SHAPELET_PHOTOM_BAD_MOMENTS;
        return _status;
    }

    Definition definition;
    definition.frames.insert(definition::Frame::make(0, *exp, fp, _bitmask, _usePixelWeights));

    _fp = definition.frames[0].getFootprint();
    if (_fp->getArea() == 0) {
        _status |= algorithms::Flags::PHOTOM_NO_FOOTPRINT;
        return _status;
    }

    //fit both a point source and galaxy model for the same object
    definition::Object galaxy = definition::Object::makeGalaxy(
        GALAXY_ID, _basis, *ellipse, 
        _isEllipticityActive, 
        _isRadiusActive, 
        _isPositionActive
    );
    definition::Object star(STAR_ID);
    //link the position object of the two models
    star.getPosition() = galaxy.getPosition();

    definition.objects.insert(galaxy);
    definition.objects.insert(star);

    _evaluator = Evaluator::make(definition);
    BruteForceSourceOptimizer optimizer;
    bool success = optimizer.solve(_evaluator, _nTestPoints);

    ndarray::Array<const double, 1, 1> coeffView;
    double fluxVar;
    ndarray::Array<const double, 2, 1> covar;
    if(success) {
        _param = ndarray::copy(optimizer.getBestParameters());        
        _coeff = ndarray::copy(optimizer.getBestCoefficients());
        covar = optimizer.getCoefficientCovariance();
        if (!optimizer.isBestSafe()) {
            _status |= algorithms::Flags::SHAPELET_PHOTOM_INVERSION_UNSAFE;
        }
        grid::Object const & object = grid::find(
            _evaluator->getGrid()->objects, GALAXY_ID
        );
        grid::Source const & source = object.sources[0];
        _ellipse.reset(new Ellipse(object.makeEllipse(_param.begin())));
        _flux = source.computeFluxMean(_param, _coeff);
        fluxVar = source.computeFluxVariance(_param, covar);

        coeffView = _coeff;
    }
    else if (_fitDeltaFunction) { 
        _status |= algorithms::Flags::SHAPELET_PHOTOM_GALAXY_FAIL;
        _status |= GALAXY_MODEL_FAILED;        
        //try fitting just the point source model
        definition.objects.erase(GALAXY_ID);
        _evaluator = Evaluator::make(definition);
        Evaluation evaluation(_evaluator);
        _coeff = ndarray::allocate(getBasisSize() + 1);
        try{
            _param = ndarray::copy(evaluation.getParameters());            
            _coeff[ndarray::view(getBasisSize(), getBasisSize() + 1)] = evaluation.getCoefficients();
            covar = evaluation.getCoefficientFisherMatrix();
        }
        catch (...){
            _status |= algorithms::Flags::SHAPELET_PHOTOM_INVERSION_FAIL;
            _ellipse.reset();
        }

        grid::Object const & object = find(
            _evaluator->getGrid()->objects, STAR_ID
        );
        _ellipse.reset( 
            new Ellipse(EllipseCore(0,0,0), object.makePoint(_param.begin()))
        ); 

        //fill basis coeffiecient with NaN
        _coeff[ndarray::view(0, getBasisSize())] = std::numeric_limits<double>::quiet_NaN();
        coeffView = _coeff[ndarray::view(getBasisSize(), getBasisSize() + 1)];
    }
    if(_fitDeltaFunction) {  
        grid::Object const & object = find(
            _evaluator->getGrid()->objects, STAR_ID
        );
        grid::Source const & source = object.sources[0];

        _flux += source.computeFluxMean(_param, coeffView);
        fluxVar += source.computeFluxVariance(_param, covar);        
    }
    _fluxErr = sqrt(fluxVar);
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
