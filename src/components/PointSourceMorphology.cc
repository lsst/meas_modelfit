// -*- lsst-c++ -*-
/**
 * @file
 * Implementation of PointSourceMorphology
 */
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"

namespace components = lsst::meas::multifit::components;

components::Morphology::Ptr components::PointSourceMorphology::create(
    boost::shared_ptr<ParameterVector const> const & linearParameters,
    ParameterConstIterator morphologyParameterIter
) const {
    return Morphology::Ptr(new PointSourceMorphology(linearParameters,morphologyParameterIter));
}

components::MorphologyProjection::Ptr components::PointSourceMorphology::makeProjection(
    lsst::afw::geom::Extent2I const & kernelDimensions,
    lsst::afw::geom::AffineTransform::ConstPtr const & transform
) const {
    return boost::make_shared<PointSourceMorphologyProjection>(
        boost::static_pointer_cast<PointSourceMorphology const>(shared_from_this()),
        kernelDimensions, transform
    );
}
