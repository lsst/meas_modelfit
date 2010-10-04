// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
 * @file
 * Support for Models which are defined by an astrometry and a morphology
 * component
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H

#include <boost/make_shared.hpp>

#include "lsst/afw/geom/Box.h"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModelProjection;

/**
 *  Derived Model defined by its a Astrometry and a Morphology components.
 *
 *  The number of nonlinear parameters is the sum of the Astrometry parameters
 *  and the Morphology parameters. When calling setNonlinearParameters on a
 *  ComponentModel, the input ParameterVector must be organized as 
 *  @sa Model
 *  @sa ComponentModelFactory
 *  @sa ComponentModelProjection
 */
class ComponentModel : public Model {
public:
    typedef boost::shared_ptr<ComponentModel> Ptr;
    typedef boost::shared_ptr<ComponentModel const> ConstPtr;

    virtual lsst::afw::detection::Footprint::Ptr computeProjectionFootprint(
        lsst::afw::detection::Psf::ConstPtr const & psf,
        lsst::afw::geom::AffineTransform const & transform
    ) const;

    virtual lsst::afw::geom::BoxD computeProjectionEnvelope(
        lsst::afw::detection::Psf::ConstPtr const & psf,
        lsst::afw::geom::AffineTransform const & transform
    ) const;

    virtual lsst::afw::geom::ellipses::Ellipse::Ptr computeBoundingEllipse() const;
    virtual lsst::afw::geom::Point2D getPosition() const {
        return _astrometry->computePosition();
    }

    virtual Model::Ptr clone() const { 
        return Model::Ptr(new ComponentModel(*this)); 
    }

    /**
     *  Immutable access to the Astrometry component.
     */    
    CONST_PTR(lsst::meas::multifit::components::Astrometry) getAstrometry() const {
        return _astrometry; 
    }

    /**
     *  Immutable access to the the Morphology component.
     */
    CONST_PTR(lsst::meas::multifit::components::Morphology) getMorphology() const { 
        return _morphology; 
    }
    
    static ComponentModel::Ptr create(
        CONST_PTR(lsst::meas::multifit::components::Astrometry) const & astrometry, 
        CONST_PTR(lsst::meas::multifit::components::Morphology) const & morphology
    ) {
        return ComponentModel::Ptr(
            new ComponentModel(astrometry, morphology, true)
        );
    }
protected:
    virtual boost::shared_ptr<ModelProjection> makeProjection(
        lsst::afw::detection::Psf::ConstPtr const & psf,
        lsst::afw::geom::AffineTransform const & transform,
        CONST_PTR(lsst::afw::detection::Footprint) const & footprint
    ) const;
private:
    explicit ComponentModel(
        lsst::meas::multifit::components::Astrometry::ConstPtr const & astrometry,
        lsst::meas::multifit::components::Morphology::ConstPtr const & morphology,
        bool initializeParameters=false
    );
    explicit ComponentModel(ComponentModel const & model);

    void _initializeFromComponents(
        lsst::meas::multifit::components::Astrometry::ConstPtr const & astrometryTemplate,
        lsst::meas::multifit::components::Morphology::ConstPtr const & morphologyTemplate,
        bool initializeParameters=false
    );

    /**
     * Iterator at the start of the nonlinear parameters relevant to the
     * astrometry component
     */
    ParameterConstIterator _beginAstrometry() const {
        return getNonlinearParameterIter();
    }
    /**
     * Iterator at the start of the nonlinear parameters relevant to the
     * morphology component
     */
    ParameterConstIterator _beginMorphology() const {
        return getNonlinearParameterIter() + _astrometry->getParameterSize();
    }


    virtual void _handleLinearParameterChange();

    virtual void _handleNonlinearParameterChange();

    lsst::meas::multifit::components::Astrometry::Ptr _astrometry;
    lsst::meas::multifit::components::Morphology::Ptr _morphology;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H
