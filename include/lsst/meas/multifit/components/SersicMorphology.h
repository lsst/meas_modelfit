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
 
#ifndef LSST_MEAS_MULTIFIT_SERSIC_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_SERSIC_MORPHOLOGY_H

#include "lsst/meas/multifit/components/Morphology.h"
#include "lsst/meas/multifit/Cache.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class SersicMorphologyProjection;

/**
 * Derived Morphology component for fitting small galaxies
 */
class SersicMorphology : public Morphology {
public:
    enum LinearParameters {
        FLUX=0,
        LINEAR_SIZE
    };
    enum NonlinearParameters {
        GAMMA1=0,
        GAMMA2,
        KAPPA,
        SERSIC_INDEX,
        NONLINEAR_SIZE
    };

    typedef boost::shared_ptr<SersicMorphology> Ptr;
    typedef boost::shared_ptr<SersicMorphology const> ConstPtr;

  
    /**
     * Named SersicMorphology constructor
     */
    static SersicMorphology::Ptr create(
        Parameter const & flux,
        lsst::afw::geom::ellipses::BaseCore const & ellipse, 
        Parameter const & sersicIndex
    ) { 
        lsst::afw::geom::ellipses::LogShear logShear(ellipse);
        boost::shared_ptr<ParameterVector> linear(new ParameterVector(LINEAR_SIZE));
        boost::shared_ptr<ParameterVector> nonlinear(new ParameterVector(NONLINEAR_SIZE));
        *linear << flux;
        *nonlinear << logShear.getVector(), sersicIndex;

        return SersicMorphology::Ptr(new SersicMorphology(linear,nonlinear));
        
    }

    Parameter const & getFlux() const {
        return *(_linearParameters->data() + FLUX);
    }
    Parameter const & getSersicIndex() const {
        return *(beginNonlinear() + SERSIC_INDEX);
    }

    virtual PTR(lsst::afw::geom::ellipses::BaseCore) computeBoundingEllipseCore() const;

    virtual MorphologyProjection::Ptr makeProjection(
        lsst::afw::geom::Extent2I const & kernelSize,
        lsst::afw::geom::AffineTransform const & transform
    ) const;

    virtual int const getNonlinearParameterSize() const {return NONLINEAR_SIZE;}

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start=0
    ) const;

    static void setSersicCache(Cache::ConstPtr const & cache) {
        _cache=cache;
    }
protected:
    friend class SersicMorphologyProjection;

    virtual void _handleNonlinearParameterChange();

    multifit::Cache::Interpolator::ConstPtr const & getInterpolator() const {
        return _interpolator;
    }
    multifit::Cache::Interpolator::ConstPtr const & getDerivativeInterpolator() const {
        return _derivativeInterpolator;
    }

    multifit::Cache::ConstPtr const & getSersicCache() const {
        return _cache;
    }

    /**
     * Construct a SersicMorphology
     *
     * For use within ComponentModel
     * @sa Morphology::create()
     */
    SersicMorphology(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start =0
    ) : Morphology(linearParameters, nonlinearParameters, start) {        
        _handleNonlinearParameterChange();
    }



private:
    Cache::Interpolator::ConstPtr _interpolator;
    Cache::Interpolator::ConstPtr _derivativeInterpolator;
    static Cache::ConstPtr _cache;
};


}}}} //namespace lsst::meas::multifit
#endif
