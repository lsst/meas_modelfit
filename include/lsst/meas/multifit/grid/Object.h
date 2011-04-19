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

namespace lsst { namespace meas { namespace multifit { namespace grid {

class Object : public definition::Object {
public:

    typedef Array<Source> SourceArray;

    Object(definition::Object const & definition_, int offset, int frameCount, int filterCount);

    Object(Object const & other) :
        definition::Object(other),
        coefficientOffset(other.coefficientOffset),
        coefficientCount(other.coefficientCount),
        sources(),
        extra(0)
    {}

    int coefficientOffset;
    int coefficientCount;
#ifndef SWIG
    SourceArray sources;
#endif
    mutable void * extra;

    PositionComponent & getPosition() const {
        return static_cast<PositionComponent &>(*this->position);
    }

    RadiusComponent & getRadius() const {
        return static_cast<RadiusComponent &>(*this->radius);
    }

    EllipticityComponent & getEllipticity() const {
        return static_cast<EllipticityComponent &>(*this->ellipticity);
    }

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

};

inline afw::geom::Point2D const Source::getReferencePoint() const {
    return transform(object.position->getReference());
}

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_Object
