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

#ifndef LSST_MEAS_MULTIFIT_ShapeletModelBasis
#define LSST_MEAS_MULTIFIT_ShapeletModelBasis

#include "lsst/meas/multifit/ModelBasis.h"
#include "lsst/afw/math/shapelets.h"
#include <boost/serialization.export.hpp>

namespace lsst { namespace meas { namespace multifit {
    class ShapeletModelBasis;
}}}

namespace boost { namespace serialization {
    template <class Archive>
    void load_construct_data(
        Archive &, lsst::meas::multifit::ShapeletModelBasis *, unsigned int const
    );
    
    template <class Archive>
    void save_construct_data(
        Archive &, const lsst::meas::multifit::ShapeletModelBasis *, unsigned int const
    );
}}
namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief An ModelBasis subclass for a single, scaled shapelet expansion.
 *
 *  ShapeletModelBasis subclasses should be immutable.
 */
class ShapeletModelBasis : public ModelBasis {
public:

    typedef boost::shared_ptr<ShapeletModelBasis> Ptr;

    static void setPsfShapeletOrder(int order) { getPsfShapeletOrderRef() = order; }
    static int getPsfShapeletOrder() { return getPsfShapeletOrderRef(); }

    /**
     *  @brief Convolve the basis with the given local PSF, returning a new basis with the same
     *         parametrization.
     */
    virtual ModelBasis::Ptr convolve(CONST_PTR(LocalPsf) const & psf) const;

    /**
     *  @brief Convolve the basis with the given ShapeletFunction, returning a new basis with the same
     *         parametrization.
     */
    ModelBasis::Ptr convolve(lsst::afw::math::shapelets::ShapeletFunction const & psf) const;

    /**
     *  @brief Convolve the basis with the given MultiShapeletFunction, returning a new basis with the same
     *         parametrization.
     */
    ModelBasis::Ptr convolve(lsst::afw::math::shapelets::MultiShapeletFunction const & psf) const;
    
    /// @brief Order of the shapelet expansion.
    int getOrder() const { return _order; };

    /// @brief Order of the shapelet expansion.
    double getScale() const { return _scale; };

    static Ptr make(int order, double scale=1.0) {
        return Ptr(new ShapeletModelBasis(order, scale));
    }

protected:

    virtual void _evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    virtual void _evaluateRadialProfile(
        lsst::ndarray::Array<Pixel,2,1> const & profile,
        lsst::ndarray::Array<Pixel const,1,1> const & radii
    ) const;

private:
    ShapeletModelBasis(int order, double scale);

#ifndef SWIG
    template <class Archive>
    friend void boost::serialization::load_construct_data(
        Archive & ar, ShapeletModelBasis * basis, unsigned const int version
    );
    template <class Archive>
    friend void boost::serialization::save_construct_data(
        Archive & ar, const ShapeletModelBasis * basis, unsigned const int version
    );

    friend class boost::serialization::access;
    template <typename Archive> void serialize(
        Archive & ar, unsigned int const version
    ){
        boost::serialization::base_object<ModelBasis>(*this);
    }
#endif

    static int & getPsfShapeletOrderRef();

    int _order;
    double _scale;
};

}}} // namespace lsst::meas::multifit

#ifndef SWIG

BOOST_CLASS_EXPORT_GUID(lsst::meas::multifit::ShapeletModelBasis, "lsst::meas::multifit::ShapeletModelBasis");

namespace boost { namespace serialization {
    
template <class Archive>
void save_construct_data(
    Archive & ar, 
    const lsst::meas::multifit::ShapeletModelBasis * basis, 
    unsigned int const version
) {
    int order = basis->getOrder();
    double scale = basis->getScale();
    ar << make_nvp("order", order);
    ar << make_nvp("scale", scale);
}
template <class Archive>
void load_construct_data(
    Archive & ar, 
    lsst::meas::multifit::ShapeletModelBasis * basis, 
    unsigned int const version
) {
    int order;
    double scale;
    ar >> make_nvp("order", order);
    ar >> make_nvp("scale", scale);
    ::new(basis)lsst::meas::multifit::ShapeletModelBasis(order, scale);
}
}}
#endif //!SWIG

#endif // !LSST_MEAS_MULTIFIT_ShapeletModelBasis
