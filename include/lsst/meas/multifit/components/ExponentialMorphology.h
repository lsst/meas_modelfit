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
 
#ifndef LSST_MEAS_MULTIFIT_EXPONENTIAL_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_EXPONENTIAL_MORPHOLOGY_H

#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class ExponentialMorphologyProjection;

/**
 * Derived Morphology component for fitting small galaxies
 */
class ExponentialMorphology : public Morphology {
public:
    enum LinearParameters {
        FLUX=0,
        LINEAR_SIZE
    };
    enum NonlinearParameters {
        GAMMA1=0,
        GAMMA2,
        KAPPA,
        NONLINEAR_SIZE
    };

    typedef boost::shared_ptr<ExponentialMorphology> Ptr;
    typedef boost::shared_ptr<ExponentialMorphology const> ConstPtr;

  
    /**
     * Named ExponentialMorphology constructor
     */
    static ExponentialMorphology::Ptr create(
        Parameter const & flux,
        lsst::afw::geom::ellipses::BaseCore const & ellipse
    ) { 
        lsst::afw::geom::ellipses::LogShear logShear(ellipse);
        boost::shared_ptr<ParameterVector> linear(new ParameterVector(LINEAR_SIZE));
        boost::shared_ptr<ParameterVector> nonlinear(new ParameterVector(NONLINEAR_SIZE));
        *linear << flux;
        *nonlinear << logShear.getVector();
        return ExponentialMorphology::Ptr(new ExponentialMorphology(linear,nonlinear));
        
    }

    Parameter const & getFlux() const {
        return *(_linearParameters->data() + FLUX);
    }

    virtual PTR(lsst::afw::geom::ellipses::BaseCore) computeBoundingEllipseCore() const;

    virtual MorphologyProjection::Ptr makeProjection(
        lsst::afw::geom::Extent2I const & kernelSize,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const;

    virtual int const getNonlinearParameterSize() const {return NONLINEAR_SIZE;}

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start=0
    ) const;

protected:
    friend class ExponentialMorphologyProjection;

    virtual void _handleNonlinearParameterChange();

    /**
     * Construct a ExponentialMorphology
     *
     * For use within ComponentModel
     * @sa Morphology::create()
     */
    ExponentialMorphology(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start =0
    ) : Morphology(linearParameters, nonlinearParameters, start) {
        _handleNonlinearParameterChange();
    }

};


}}}} //namespace lsst::meas::multifit
#endif
