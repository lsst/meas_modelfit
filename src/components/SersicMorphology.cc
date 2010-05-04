#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
#include "lsst/afw/geom/ellipses/LogShear.h"

namespace components = lsst::meas::multifit::components;

lsst::afw::geom::ellipses::Core::Ptr 
components::SersicMorphology::computeBoundingEllipseCore() const {  
    ParameterConstIterator params(beginNonlinear());
    return boost::make_shared<lsst::afw::geom::ellipses::LogShear> (
        params[GAMMA1], params[GAMMA2], params[KAPPA]
    );
}
components::Morphology::Ptr components::SersicMorphology::create(
    boost::shared_ptr<ParameterVector const> const & linearParameters,
    boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
    size_t const & start
) const {
    return Morphology::Ptr(
        new SersicMorphology(linearParameters, nonlinearParameters, start)
    );
}

components::MorphologyProjection::Ptr components::SersicMorphology::makeProjection(
    lsst::afw::geom::Extent2I const & kernelDimensions,
    lsst::afw::geom::AffineTransform::ConstPtr const & transform
) const {
    return SersicMorphologyProjection::Ptr(new SersicMorphologyProjection(
        boost::static_pointer_cast<SersicMorphology const>(shared_from_this()),
        kernelDimensions,
        transform
    ));
}

void components::SersicMorphology::_handleNonlinearParameterChange() {
    _indexFunctor = SersicCache::getInstance()->getRowFunctor(getSersicIndex());
}
