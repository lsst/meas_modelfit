#ifndef LSST_MEAS_MULTIFIT_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_MODEL_FACTORY_H

#include "lsst/afw/geom/Point.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"

namespace lsst{
namespace meas{
namespace multifit {

class ModelFactory {
public:
    static Model::Ptr createPointSourceModel(
        Parameter const & flux,
        lsst::afw::geom::Point2D const & centroid
    ) {
        components::Morphology::Ptr morphology = 
            components::PointSourceMorphology::create(flux);
        components::Astrometry::Ptr astrometry(
            new components::Astrometry(centroid)
        );
        return ComponentModel::create(astrometry, morphology);
    }

    static Model::Ptr createSersicModel(
        Parameter const & flux,
        lsst::afw::geom::Point2D const & centroid,
        lsst::afw::geom::ellipses::Core const & ellipse,
        Parameter const & sersicIndex
    ) {
        components::Morphology::Ptr morphology = 
            components::SersicMorphology::create(flux, ellipse, sersicIndex);
        components::Astrometry::Ptr astrometry(
            new components::Astrometry(centroid)
        );
        return ComponentModel::create(astrometry, morphology); 
    }
};

}}}

#endif
