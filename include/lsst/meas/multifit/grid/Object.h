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

    SourceArray sources;

    /// @brief The number of coefficients for this object per Frame.
    int const getSourceCoefficientCount() const {
        return getBasis() ? getBasis()->getSize() : 1;
    }

    /// @brief The offset of this object's coefficients in the grid's full coefficient array.
    int const getCoefficientOffset() const { return _coefficientOffset; }

    /// @brief The total number of coefficients for this object.
    int const getCoefficientCount() const { return _coefficientCount; }

#ifndef SWIG
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

    /// @brief Throw an exception (LogicErrorException) if the object lacks a radius or ellipticity.
    void requireEllipse() const;

    /// @brief Construct the point corresponding to this object from a parameter vector.
    lsst::afw::geom::Point2D makePoint(double const * paramIter) const;

    /**
     *  @brief Fill the elements of a parameter vector with values from the given point.
     *
     *  If the position component of the object is inactive, no action is taken.
     */
    void readPoint(double * paramIter, lsst::afw::geom::Point2D const & point) const;

    /**
     *  @brief Given a symmetric matrix corresponding to a full parameter vector, extract
     *         a 2x2 matrix corresponding to the point.
     *
     *  The returned matrix will be zero if the position is inactive.
     *
     *  The ordering of ellipse parameters is (e1, e2, r, x, y).
     */
    Eigen::Matrix2d extractPointMatrix(Eigen::MatrixXd const & matrix) const;

    /**
     *  @brief Given a symmetric matrix corresponding to the full parameter vector, set the
     *         2x2 block corresponding to the point.
     *
     *  If the position is inactive, the full matrix will not be modified.
     *
     *  The ordering of ellipse parameters is (e1, e2, r, x, y).
     */
     void insertPointMatrix(Eigen::MatrixXd & full, Eigen::Matrix2d const & block) const;

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

    /// @brief Fill the elements of a parameter vector with values from the given ellipse.
    void readEllipse(double * paramIter, lsst::afw::geom::ellipses::Ellipse const & ellipse) const;

    /**
     *  @brief Given a symmetric matrix corresponding to a full parameter vector, extract
     *         a 5x5 matrix corresponding to the ellipse.
     *
     *  Rows and columns that correspond to inactive parameters will be zero.
     */
    Eigen::Matrix5d extractEllipseMatrix(Eigen::MatrixXd const & matrix) const;

    /**
     *  @brief Given a symmetric matrix corresponding to the full parameter vector, set the
     *         5x5 block corresponding to the ellipse.
     *
     *  Block rows and columns that correspond to inactive parameters will be ignored.
     */
    void insertEllipseMatrix(Eigen::MatrixXd & full, Eigen::Matrix5d const & block) const;

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

#ifndef SWIG
std::ostream & operator<<(std::ostream & os, Object const & obj);

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
#endif

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_Object
