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

#ifndef LSST_MEAS_MULTIFIT_SimpleInterpreter
#define LSST_MEAS_MULTIFIT_SimpleInterpreter

#include "lsst/meas/multifit/BaseInterpreter.h"
#include "lsst/meas/multifit/SimpleDistribution.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for interpreters based on SimpleDistribution.
 */
class SimpleInterpreter : public virtual BaseInterpreter {
public:

    typedef boost::shared_ptr<SimpleInterpreter> Ptr;
    typedef boost::shared_ptr<SimpleInterpreter const> ConstPtr;

    SimpleDistribution::Ptr getTarget() {
        return boost::static_pointer_cast<SimpleDistribution>(_getTarget());
    }

    SimpleDistribution::ConstPtr getTarget() const {
        return boost::static_pointer_cast<SimpleDistribution const>(_getTarget());
    }

    ///@{
    /// @brief BaseInterpreter interface implementations.
    virtual lsst::afw::geom::Point2D computePointMean(ID id) const { return extractPointMu(id); }
    virtual Eigen::Matrix2d computePointCovariance(ID id) const { return extractPointSigma(id); }
    virtual Ellipse computeEllipseMean(ID id) const { return extractEllipseMu(id); }
    virtual Eigen::Matrix5d computeEllipseCovariance(ID id) const { return extractEllipseSigma(id); }
    ///@}

    /**
     *  @brief Extract the position corresponding to the mu vector section for the given object.
     *
     *  @note while the position-related elements of the full mu vector are defined as offsets
     *        relative to the reference points of the PositionComponent objects in the grid,
     *        the points returned here are absolute positions.
     */
    lsst::afw::geom::Point2D extractPointMu(ID id) const;

    /**
     *  @brief Set the position corresponding to the mu vector section for the given object.
     *
     *  @sa extractPointMu
     */
    void insertPointMu(ID id, lsst::afw::geom::Point2D const & mu);

    /**
     *  @brief Extract the 2x2 block of the sigma matrix corresponding to the position of the given object.
     *
     *  If the position of the given object is inactive the returned matrix will be the zero matrix.
     */
    Eigen::Matrix2d extractPointSigma(ID id) const;

    /**
     *  @brief Set the 2x2 block of the sigma matrix corresponding to the position of the given object.
     *
     *  If the position of the given object is inactive no action is taken.
     */
    void insertPointSigma(ID id, Eigen::Matrix2d const & sigma);

    /// @brief Extract the ellipse corresponding to mu vector section for the given object.
    Ellipse extractEllipseMu(ID id) const;

    /**
     *  @brief Set the ellipse corresponding to the mu vector section for the given object.
     *
     *  The ellipse core will be converted to multifit::EllipseCore before assigning to the vector.
     *  Inactive parameters will be ignored.
     */
    void insertEllipseMu(ID id, Ellipse const & mu);

    /**
     *  @brief Extract the 5x5 block of the sigma matrix corresponding to the ellipse of the given object.
     *
     *  The ordering and definition of the matrix elements corresponds to the elements of 
     *  Ellipse::getParameterVector(), with core set by multifit::EllipseCore: (e1, e2, r, x, y).
     *
     *  Ellipse components that are not active will have zero variance, resulting in a singular matrix.
     */
    Eigen::Matrix5d extractEllipseSigma(ID id) const;

    /**
     *  @brief Set the 5x5 block of the sigma matrix corresponding to the ellipse of the given object.
     *
     *  Rows and columns that correspond to inactive Ellipse components will be ignored.
     *
     *  @sa extractEllipseSigma.
     */
    void insertEllipseSigma(ID id, Eigen::Matrix5d const & sigma);

protected:

    explicit SimpleInterpreter(
        Grid::Ptr const & grid, SimpleDistribution::Ptr const & target
    ) :
        BaseInterpreter(grid), _target(target) {}

    explicit SimpleInterpreter(
        Grid::Ptr const & grid, 
        SimpleDistribution::ConstPtr const & target
    ) :
        BaseInterpreter(grid), _target(boost::const_pointer_cast<SimpleDistribution>(target)) {}

    virtual BaseDistribution::Ptr _getTarget() { return _target; }
    virtual BaseDistribution::ConstPtr _getTarget() const { return _target; }

    Eigen::VectorXd & getMuRef() { return _target->_mu; }
    Eigen::MatrixXd & getSigmaRef() { return _target->_sigma; }

    Eigen::VectorXd const & getMuCRef() const { return _target->_mu; }
    Eigen::MatrixXd const & getSigmaCRef() const { return _target->_sigma; }

    void invalidateTarget() { _target->invalidate(); }

    SimpleDistribution::Ptr _target;
};

/**
 *  @brief NestedInterpreter for SimpleDistribution.
 */
class NestedSimpleInterpreter : public SimpleInterpreter, public NestedInterpreter {
public:

    typedef boost::shared_ptr<NestedSimpleInterpreter> Ptr;
    typedef boost::shared_ptr<NestedSimpleInterpreter const> ConstPtr;

    /// @brief Construct a mutable interpreter.
    static Ptr make(SimpleDistribution::Ptr const & target, Grid::Ptr const & grid) {
        return Ptr(new NestedSimpleInterpreter(grid, target));
    }

    /// @brief Construct a const interpreter.
    static ConstPtr make(SimpleDistribution::ConstPtr const & target, Grid::Ptr const & grid) {
        return ConstPtr(new NestedSimpleInterpreter(grid, target));
    }

private:

    /// @brief Construct from a multifit grid and mutable SimpleDistribution pointer.
    NestedSimpleInterpreter(
        Grid::Ptr const & grid,
        SimpleDistribution::Ptr const & target
    ) : BaseInterpreter(grid), SimpleInterpreter(grid, target), NestedInterpreter(grid) {
        ensureCompatibility();
    }

    /// @brief Construct from a multifit grid and const SimpleDistribution pointer.
    NestedSimpleInterpreter(
        Grid::Ptr const & grid,
        SimpleDistribution::ConstPtr const & target
    ) : BaseInterpreter(grid), SimpleInterpreter(grid, target), NestedInterpreter(grid) {
        ensureCompatibility();
    }

    void ensureCompatibility();

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SimpleInterpreter
