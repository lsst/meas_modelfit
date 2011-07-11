#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

SourceComponent::SourceComponent(
    Frame const & frame_, ObjectComponent const & object_, 
    CONST_PTR(afw::image::Wcs) const & wcs
) :
    frame(frame_), object(object_), 
    _transform()
{
    if (wcs) {
        if (!frame.getWcs()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "If the definition WCS is set, all frames must have a WCS."
            );
        }
	afw::geom::Point2D point = object.getPosition()->getValue();
        _transform = frame.getWcs()->linearizeSkyToPixel(point) * wcs->linearizePixelToSky(point);
    } else {
        if (frame.getWcs()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "If the definition WCS is not set, inividual frames may not have a WCS."
            );
        }
    }
    if (frame.getPsf()) {
        _localPsf = frame.getPsf()->getLocalPsf(_transform(object.getPosition()->getValue()));
    }
    if (object.getBasis()) {
        if (frame.getPsf()) {            
            _basis = object.getBasis()->convolve(_localPsf);
        } else {
            _basis = object.getBasis();
        }
    } else {
        if (!frame.getPsf()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "All objects must have a basis if any frames do not have a PSF."
            );
        }
    }
}

void SourceComponent::fillIntegration(lsst::ndarray::Array<Pixel,1,1> const & integration) const {
    if (object.getBasis()) {
        integration[ndarray::view(getCoefficientOffset(), getCoefficientOffset() + getCoefficientCount())]
            = object.getBasis()->getIntegration();
    } else {
        integration[getCoefficientOffset()] = 1.0;
    }
}

double SourceComponent::computeFlux(
    lsst::ndarray::Array<Pixel const,1,1> const & integration,
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients
) {
    return ndarray::viewAsEigen(integration).dot(ndarray::viewAsEigen(coefficients));
}

double SourceComponent::computeFluxVariance(
    lsst::ndarray::Array<Pixel const,1,1> const & integration,
    lsst::ndarray::Array<Pixel const,2,1> const & covariance
) {
    return ndarray::viewAsEigen(integration).dot(
        ndarray::viewAsEigen(covariance) * ndarray::viewAsEigen(integration)
    );
}

double SourceComponent::computeFlux(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients
) const {
    if (object.getBasis()) {
        return ndarray::viewAsEigen(
            coefficients[
                ndarray::view(getCoefficientOffset(), getCoefficientOffset() + getCoefficientCount())
            ]).dot(ndarray::viewAsEigen(object.getBasis()->getIntegration()));
    }
    return coefficients[getCoefficientOffset()];
}

double SourceComponent::computeFluxVariance(
    lsst::ndarray::Array<Pixel const,2,1> const & covariance
) const {
    if (object.getBasis()) {
        Eigen::MatrixXd sigma = ndarray::viewAsEigen(covariance).block(
            getCoefficientOffset(),
            getCoefficientOffset(),
            getCoefficientCount(),
            getCoefficientCount()
        );
        ndarray::Array<Pixel const,1,1> integration = object.getBasis()->getIntegration();
        return ndarray::viewAsEigen(integration).dot(sigma * ndarray::viewAsEigen(integration));
    }
    return covariance[getCoefficientOffset()][getCoefficientOffset()];
}

}}}} // namespace lsst::meas::multifit::grid
