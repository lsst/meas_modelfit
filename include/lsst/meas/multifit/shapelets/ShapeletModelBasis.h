// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletModelBasis
#define LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletModelBasis

#include "lsst/meas/multifit/ModelBasis.h"

namespace lsst { namespace meas { namespace multifit { namespace shapelets {

class ConvolvedShapeletModelBasis;

/**
 *  @brief An ModelBasis subclass for a single, scaled shapelet expansion.
 *
 *  ShapeletModelBasis subclasses should be immutable.
 */
class ShapeletModelBasis : public ModelBasis {
public:

    typedef boost::shared_ptr<ShapeletModelBasis> Ptr;
    typedef ConvolvedShapeletModelBasis ConvolvedBasis;

    /**
     *  @brief Convolve the basis with the given local PSF, returning a new basis with the same
     *         parametrization.
     */
    PTR(ConvolvedBasis) convolve(PTR(LocalPsf) const & psf) const {
        return boost::static_pointer_cast<ConvolvedBasis>(_convolve(psf));
    }

    /// @brief Order of the shapelet expansion.
    int getOrder() const { return _order; };

    /// @brief Order of the shapelet expansion.
    double getScale() const { return _scale; };

    static Ptr make(int order, double scale=1.0) {
        return boost::make_shared<ShapeletModelBasis>(order, scale);
    }

protected:

    virtual void _evaluate(
        lsst::ndarray::Array<double, 2, 1> const & matrix,
        PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    virtual PTR(ModelBasis) _convolve(PTR(LocalPsf) const & psf) const;

private:

    FRIEND_MAKE_SHARED_2(ShapeletModelBasis, int, double);

    ShapeletModelBasis(int order, double scale);

    int _order;
    double _scale;

    void operator=(ShapeletModelBasis const &) {}
};

}}}} // namespace lsst::meas::multifit::shapelets

#endif // !LSST_MEAS_MULTIFIT_SHAPELETS_ShapeletModelBasis
