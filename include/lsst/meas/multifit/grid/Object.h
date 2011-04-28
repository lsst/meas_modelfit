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

#ifndef LSST_MEAS_MULTIFIT_GRID_Object
#define LSST_MEAS_MULTIFIT_GRID_Object

#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/multifit/definition/Object.h"
#include "lsst/meas/multifit/grid/parameters.h"
#include "lsst/meas/multifit/grid/Array.h"
#include "lsst/meas/multifit/grid/Source.h"
#include "lsst/meas/multifit/grid/Frame.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

/**
 *  @brief An immutable and expanded version of definition::Object used in a multifit Grid.
 */
class Object : public detail::ObjectBase, private boost::noncopyable {
public:

    typedef Array<Source> SourceArray;
    
    // We can show all the internals except the parameter components because users will only ever 
    // see const Objects, but we can't let them have non-const pointers to the parameter components.
    // We'd also like to cast the components to their grid versions.

    SourceArray sources;

    /// @brief The number of coefficients for this object per Frame.
    int const getSourceCoefficientCount() const {
        return getBasis() ? getBasis()->getSize() : 1;
    }

    /// @brief The offset of this object's coefficients in the grid's full coefficient array.
    int const getCoefficientOffset() const { return _coefficientOffset; }

    /// @brief The total number of coefficients for this object.
    int const getCoefficientCount() const { return _coefficientCount; }

    //@{
    /// @brief Return the parameter components.
    PositionComponent::Ptr const getPosition() const { return _position; }
    RadiusComponent::Ptr const getRadius() const { return _radius; }
    EllipticityComponent::Ptr const getEllipticity() const { return _ellipticity; }

    template <ParameterType E>
    typename ParameterComponent<E>::Ptr const getComponent() const {
        typename ParameterComponent<E>::Ptr const * p;
        getComponentImpl(p);
        return *p;
    }
    //@}

#ifndef SWIG
    /// @brief Construct the point corresponding to this object from a parameter vector.
    lsst::afw::geom::Point2D makePoint(double const * paramIter) const;

    /**
     *  Perturb a point by changing the nth position parameter.
     *  Returns the offset of the nth position parameter and the
     *  size of the perturbation.
     */
    std::pair<int,double> perturbPoint(lsst::afw::geom::Point2D & point, int n) const;

    /**
     *  Unperturb a point by changing the nth position parameter.
     */
    void unperturbPoint(lsst::afw::geom::Point2D & point, int n, double perturbation) const {
        point[n] -= perturbation;
    }

    /// @brief Construct the ellipse corresponding to this object from a parameter vector.
    lsst::afw::geom::ellipses::Ellipse makeEllipse(double const * paramIter) const;

    /**
     *  Perturb an ellipse by changing the nth ellipse parameter.
     *  Returns the offset of the nth ellipse parameter and the
     *  size of the perturbation.
     */
    std::pair<int,double> perturbEllipse(lsst::afw::geom::ellipses::Ellipse & ellipse, int n) const;

    /**
     *  Unperturb a point by changing the nth position parameter.
     */
    void unperturbEllipse(lsst::afw::geom::ellipses::Ellipse & ellipse, int n, double perturbation) const;
#endif

private:

    friend class grid::Initializer;

    Object(definition::Object const & def, int coefficientOffset, int filterCount, int frameCount);

    void validate() const;

    void getComponentImpl(PositionComponent::Ptr const * & p) const { p = &_position; }
    void getComponentImpl(RadiusComponent::Ptr const * & p) const { p = &_radius; }
    void getComponentImpl(EllipticityComponent::Ptr const * & p) const { p = &_ellipticity; }

    int _coefficientOffset;
    int _coefficientCount;
    PositionComponent::Ptr _position;
    RadiusComponent::Ptr _radius;
    EllipticityComponent::Ptr _ellipticity;
};

inline afw::geom::Point2D const Source::getReferencePoint() const {
    return _transform(object.getPosition()->getValue());
}

inline int const Source::getCoefficientOffset() const {
    return object.getCoefficientOffset()
        + object.getSourceCoefficientCount()
        * (object.isVariable() ? frame.getFrameIndex() : frame.getFilterIndex());
}

inline int const Source::getCoefficientCount() const {
    return object.getSourceCoefficientCount();
}

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_Object
