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
 
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FIXED_NONLINEAR_MORPHOLOGY
#define LSST_MEAS_MULTIFIT_COMPONENTS_FIXED_NONLINEAR_MORPHOLOGY

#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst{
namespace meas {
namespace multifit {
namespace components {

class FixedNonlinearMorphology : public Morphology {
public:
    typedef boost::shared_ptr<FixedNonlinearMorphology> Ptr;
    typedef boost::shared_ptr<FixedNonlinearMorphology const> ConstPtr;

    static Ptr create(Morphology const & base) {
        boost::shared_ptr<ParameterVector> linear(
            new ParameterVector(base.getLinearParameters())
        ); 
        boost::shared_ptr<ParameterVector> nonlinear;
        if (base.getNonlinearParameterSize() != 0) {
            nonlinear.reset(
                new ParameterVector(base.getNonlinearParameterSize())
            );
            std::copy(
                base.beginNonlinear(), 
                base.endNonlinear(), 
                nonlinear->data()
            );            
        } else {
            nonlinear.reset(new ParameterVector());
        }
        
        Ptr fixed(new FixedNonlinearMorphology(base, linear, nonlinear));
        return fixed;
    }
    virtual lsst::afw::geom::ellipses::Core::Ptr computeBoundingEllipseCore() const {
        return _base->computeBoundingEllipseCore();
    };

    virtual int const getNonlinearParameterSize() const {return 0;}
    
    virtual MorphologyProjection::Ptr makeProjection (
        lsst::afw::geom::Extent2I const & kernelSize,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const {
        return _base->makeProjection(kernelSize, transform);
    }

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const &,
        boost::shared_ptr<ParameterVector const> const &,
        size_t const &
    ) const {
        Ptr fixed(new FixedNonlinearMorphology(*_base, _linearParameters, _nonlinearParameters));
        return fixed;
    }
protected:    
    Morphology::Ptr _base;

    virtual void _handleLinearParameterChange() {_base->_handleLinearParameterChange();}
    virtual void _handleNonlinearParameterChange() {};

    FixedNonlinearMorphology(
        Morphology const & base,
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters
    ) : Morphology(linearParameters, nonlinearParameters),
        _base(base.create(_linearParameters, _nonlinearParameters)) 
    {}
};

}}}}
#endif
