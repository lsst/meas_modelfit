#include "lsst/meas/multifit/SimpleInterpreter.h"
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

lsst::afw::geom::Point2D SimpleInterpreter::extractPointMu(ID id) const {
    grid::Object const & obj = grid::find(_grid->objects, id);
    return obj.makePoint(getMuCRef().data());
}

void SimpleInterpreter::insertPointMu(ID id, lsst::afw::geom::Point2D const & mu) {
    grid::Object const & obj = grid::find(_grid->objects, id);
    obj.readPoint(getMuRef().data(), mu);
}

Eigen::Matrix2d SimpleInterpreter::extractPointSigma(ID id) const {
    grid::Object const & obj = grid::find(_grid->objects, id);
    return obj.extractPointMatrix(getSigmaCRef());
}

void SimpleInterpreter::insertPointSigma(ID id, Eigen::Matrix2d const & sigma) {
    grid::Object const & obj = grid::find(_grid->objects, id);
    invalidateTarget();
    obj.insertPointMatrix(getSigmaRef(), sigma);
}

Ellipse SimpleInterpreter::extractEllipseMu(ID id) const {
    grid::Object const & obj = grid::find(_grid->objects, id);
    return obj.makeEllipse(getMuCRef().data());
}

void SimpleInterpreter::insertEllipseMu(ID id, Ellipse const & mu) {
    grid::Object const & obj = grid::find(_grid->objects, id);
    obj.readEllipse(getMuRef().data(), mu);
}

Eigen::Matrix5d SimpleInterpreter::extractEllipseSigma(ID id) const {
    grid::Object const & obj = grid::find(_grid->objects, id);
    return obj.extractEllipseMatrix(getSigmaCRef());
}

void SimpleInterpreter::insertEllipseSigma(ID id, Eigen::Matrix5d const & sigma) {
    grid::Object const & obj = grid::find(_grid->objects, id);
    invalidateTarget();
    obj.insertEllipseMatrix(getSigmaRef(), sigma);
}

double NestedSimpleInterpreter::computeFluxMean(ID object, ID frame) const {
    throw LSST_EXCEPT(
        lsst::pex::exceptions::LogicErrorException,
        "Flux calculations not supported by NestedInterpreters."
    );
}

double NestedSimpleInterpreter::computeFluxVariance(ID object, ID frame) const {
    throw LSST_EXCEPT(
        lsst::pex::exceptions::LogicErrorException,
        "Flux calculations not supported by NestedInterpreters."
    );
}

Eigen::VectorXd NestedSimpleInterpreter::computeParameterMean() const { return getMuCRef(); }

Eigen::VectorXd NestedSimpleInterpreter::computeCoefficientMean() const {
    throw LSST_EXCEPT(
        lsst::pex::exceptions::LogicErrorException,
        "NestedSimpleInterpreter::computeCoefficientMean not yet implemented."
    );    
}

void NestedSimpleInterpreter::ensureCompatibility() {
    if (_grid->getParameterCount() != _target->getDimensionality()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Grid parameter size (%d) does not match distribution dimensionality (%d).")
             % _grid->getParameterCount() % _target->getDimensionality()).str()
        );
    }
}

Eigen::VectorXd UnifiedSimpleInterpreter::computeParameterMean() const {
    if (_grid->getParameterCount() == 0)
        return Eigen::VectorXd();
    return getMuCRef().segment(0, _grid->getParameterCount());
}

Eigen::VectorXd UnifiedSimpleInterpreter::computeCoefficientMean() const {
    return getMuCRef().segment(_grid->getParameterCount(), _grid->getCoefficientCount());
}

double UnifiedSimpleInterpreter::computeFluxMean(ID objectId, ID frameId) const {
    grid::Object const & object = grid::find(_grid->objects, objectId);
    grid::Frame const & frame = grid::find(_grid->frames, frameId);
    grid::Source const & source = object.sources[frame.getFrameIndex()];
    assert(&object == &source.object);
    assert(&frame == &source.frame);
    if (object.getBasis()) {
        Eigen::VectorXd coeff = getMuCRef().segment(
            _grid->getParameterCount() + source.getCoefficientOffset(),
            source.getCoefficientCount()
        );
        Ellipse ellipse(object.makeEllipse(getMuCRef().data()));
        ndarray::Array<Pixel,1,1> integration(ndarray::allocate(object.getBasis()->getSize()));
        object.getBasis()->integrate(integration);
        return ellipse.getCore().getArea() * ndarray::viewAsEigen(integration).dot(coeff) / M_PI;
    }
    return getMuCRef()[_grid->getParameterCount() + source.getCoefficientOffset()];
}

double UnifiedSimpleInterpreter::computeFluxVariance(ID objectId, ID frameId) const {
    grid::Object const & object = grid::find(_grid->objects, objectId);
    grid::Frame const & frame = grid::find(_grid->frames, frameId);
    grid::Source const & source = object.sources[frame.getFrameIndex()];
    assert(&object == &source.object);
    assert(&frame == &source.frame);
    if (object.getBasis()) {
        Eigen::VectorXd sigma = getSigmaCRef().block(
            _grid->getParameterCount() + source.getCoefficientOffset(),
            _grid->getParameterCount() + source.getCoefficientOffset(),
            source.getCoefficientCount(),
            source.getCoefficientCount()
        );
        Ellipse ellipse(object.makeEllipse(getMuCRef().data()));
        ndarray::Array<Pixel,1,1> integration(ndarray::allocate(object.getBasis()->getSize()));
        object.getBasis()->integrate(integration);
        // TODO: factor in ellipse covariance
        double result = ellipse.getCore().getArea() / M_PI; result *= result;
        result *= ndarray::viewAsEigen(integration).dot(sigma * ndarray::viewAsEigen(integration));
        return result;
    }
    return getSigmaCRef()(
        _grid->getParameterCount() + source.getCoefficientOffset(),
        _grid->getParameterCount() + source.getCoefficientOffset()
    );
}

void UnifiedSimpleInterpreter::ensureCompatibility() {
    int nUnified = _grid->getParameterCount() + _grid->getCoefficientCount();
    if (nUnified != _target->getDimensionality()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Grid unified parameter size (%d) does not match distribution dimensionality (%d).")
             % nUnified % _target->getDimensionality()).str()
        );
    }
}

}}} // namespace lsst::meas::multifit
