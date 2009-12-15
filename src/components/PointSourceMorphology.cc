#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"

namespace components = multifit::components;

components::Morphology::Ptr components::PointSourceMorphology::create(
    boost::shared_ptr<ParameterVector const> const & linearParameters,
    ParameterConstIterator morphologyParameterIter
) const {
    return Morphology::Ptr(new PointSourceMorphology(linearParameters,morphologyParameterIter));
}

components::MorphologyProjection::Ptr components::PointSourceMorphology::makeProjection(
    int kernelSize,
    lsst::afw::geom::AffineTransform::ConstPtr const & transform
) const {
    return boost::make_shared<PointSourceMorphologyProjection>(
        boost::static_pointer_cast<PointSourceMorphology const>(shared_from_this()),
        kernelSize, transform
    );
}
