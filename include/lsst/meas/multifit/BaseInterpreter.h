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

#ifndef LSST_MEAS_MULTIFIT_BaseInterpreter
#define LSST_MEAS_MULTIFIT_BaseInterpreter

#include "lsst/meas/multifit/Grid.h"
#include "lsst/meas/multifit/BaseDistribution.h"
 
namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for "interpreters" that relate probability distributions to the parameters
 *         and coefficients defined by a multifit grid.
 *
 *  The interpreters mostly form a classically diamond-shaped multiple-inheritance hierarchy.  One base
 *  of each interpreter sets the class of distributions the interpreter operates on, while the other
 *  sets how to relate the distribution parameters to the multifit grid.  Most of the functionality
 *  and interface of an interpreter only exists at the leaf class level.
 *
 *  The only state an interpreter has should be the multifit grid and the distribution, which should 
 *  be held by a shared_ptr.  Interpreters should only exist inside shared_ptrs, and should be held in
 *  a shared_ptr-to-const when constructed with a distribution shared_ptr-to-const.  Because they
 *  are noncopyable, this allows interpreter constness to apply to their held distribution even though
 *  the distribution is shared.
 *
 *  Many of the virtual subclasses of BaseInterpreter don't actually do anything yet; they're largely
 *  placeholders.  And because we'd like to preserve the diamond inheritance structure for all
 *  leaf-class interpeters, some of them will remain empty even if others do define an interface.
 */
class BaseInterpreter : private boost::noncopyable {
public:

    typedef boost::shared_ptr<BaseInterpreter> Ptr;
    typedef boost::shared_ptr<BaseInterpreter const> ConstPtr;

    /// @brief Return the distribution being interpreted.
    BaseDistribution::Ptr getTarget() { return _getTarget(); }

    /// @brief Return the distribution being interpreted (const).
    BaseDistribution::ConstPtr getTarget() const { return _getTarget(); }

    /// @brief Return the mean position of the object with the given ID.
    virtual lsst::afw::geom::Point2D computePointMean(ID id) const = 0;

    /// @brief Return the covariance matrix of the position of the object with the given ID.
    virtual Eigen::Matrix2d computePointCovariance(ID id) const = 0;

    /// @brief Return the mean ellipse of the object with the given ID.
    virtual Ellipse computeEllipseMean(ID id) const = 0;

    /// @brief Return the covariance matrix of the ellipse of the object with the given ID.
    virtual Eigen::Matrix5d computeEllipseCovariance(ID id) const = 0;

    virtual ~BaseInterpreter() {}

protected:
    
    virtual BaseDistribution::Ptr _getTarget() = 0;
    virtual BaseDistribution::ConstPtr _getTarget() const = 0;

};

/**
 *  @brief The distribution models the grid parameters only, with grid 
 *         coefficients marginalized out or fixed.
 */
class ParameterInterpreter : public virtual BaseInterpreter {
public:

    typedef boost::shared_ptr<ParameterInterpreter> Ptr;
    typedef boost::shared_ptr<ParameterInterpreter const> ConstPtr;

};

/**
 *  @brief An interpreter for distributions that model the concatenation of the
 *         multifit grid's parameters and coefficients, representing their joint distribution.
 */
class UnifiedInterpreter : public virtual BaseInterpreter {
public:

    typedef boost::shared_ptr<UnifiedInterpreter> Ptr;
    typedef boost::shared_ptr<UnifiedInterpreter const> ConstPtr;

    // TODO

};

/**
 *  @brief The distribution models the joint distribution of grid parameters and coefficients
 *         by nesting the conditional distribution of the coefficients within the marginal distribution
 *         of the parameters.
 */
class NestedInterpreter : public virtual BaseInterpreter {
public:

    typedef boost::shared_ptr<NestedInterpreter> Ptr;
    typedef boost::shared_ptr<NestedInterpreter const> ConstPtr;

    // TODO

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseInterpreter
