// -*- LSST-C++ -*-
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

#ifndef LSST_MEAS_MULTIFIT_SHAPELETS_ConvolvedShapeletModelBasis
#define LSST_MEAS_MULTIFIT_SHAPELETS_ConvolvedShapeletModelBasis

#include "lsst/meas/multifit/shapelets/ShapeletModelBasis.h"
#include "lsst/meas/multifit/shapelets/ShapeletConvolution.h"

namespace lsst { namespace meas { namespace multifit { namespace shapelets {

/**
 *  @brief A convolved ShapeletModelBasis.
 */
class ConvolvedShapeletModelBasis : public ModelBasis {
public:

    typedef boost::shared_ptr<ConvolvedShapeletModelBasis> Ptr;
    typedef ModelBasis ConvolvedBasis;

    int getOrder() const { return _convolution->getColOrder(); }

    double getScale() const { return _scale; }

protected:

    FRIEND_MAKE_SHARED_2(
        ConvolvedShapeletModelBasis,
        lsst::meas::multifit::shapelets::ShapeletModelBasis,
        lsst::afw::math::shapelets::EllipticalShapeletFunction
    );

    explicit ConvolvedShapeletModelBasis(
        ShapeletModelBasis const & basis,
        lsst::afw::math::shapelets::EllipticalShapeletFunction const & psf
    );

    virtual void _evaluate(
        lsst::ndarray::Array<double, 2, 1> const & matrix,
        PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    virtual ModelBasis::Ptr _convolve(PTR(LocalPsf) const & psf) const { return ModelBasis::Ptr(); }

private:

    ShapeletConvolution::Ptr _convolution;
    ShapeletModelBasis::Ptr _frontBasis;
    double _scale;
};

}}}} // namespace lsst::meas::multifit::shapelets

#endif // !LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletModelBasis
