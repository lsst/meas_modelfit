// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_parameters_h_INCLUDED
#define LSST_MEAS_MULTIFIT_parameters_h_INCLUDED

#include "lsst/base.h"
#include "lsst/afw/geom/ellipses/BaseCore.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Class that converts between two related (usually equivalent) parametrizations of a model.
 *
 *  Always obtained via ParameterDefinition::makeConversionTo().
 */
class ParameterConverter : private boost::noncopyable {
public:

    /// Convert from the input parametrization to the output parametrization
    virtual void apply(
        ndarray::Array<double const,1,1> const & input,
        ndarray::Array<double,1,1> const & output
    ) const = 0;

    /**
     *  @brief Return the determinant of the Jacobian from input to output.
     *
     *  More precisely, if @f$\theta_{\text{out}}=f(\theta_{\text{in}})@f$, then
     *  this function returns:
     *  @f[
     *    J = \left|\frac{\partial f}{\partial \theta_{\text{in}}}\right|
     *  @f]
     */
    virtual double computeJacobian(ndarray::Array<double const,1,1> const & input) const = 0;

    virtual ~ParameterConverter() {}

};

/**
 *  @brief Class that defines how to intepret a nonlinear parameter vector.
 *
 *  All ParameterDefinition classes are named singletons, stored in a global registry.
 *  The registration is performed by the base class constructor, so each derived class
 *  should instantiate itself exactly once in a static-scope variable.
 *
 *  All afw::geom::ellipse::BaseCore subclass names are valid size-3 ParameterDefinitions,
 *  and will be registered on first use.
 */
class ParameterDefinition : private boost::noncopyable {
public:

    std::string const name;   ///< Globally unique name of this parameter definition
    int const size;           ///< Number of parameters

    //@{
    /**
     *  Comparisons with other parameter definitions; equivalent to pointer comparison because all
     *  instances are singletons.
     */
    bool operator==(ParameterDefinition const & other) const { return this == &other; }
    bool operator!=(ParameterDefinition const & other) const { return !(*this == other); }
    //@}

    /**
     *  @brief Obtain an object that can convert from this definition to another.
     *
     *  If the conversion is impossible, a null pointer will be returned.
     */
    virtual PTR(ParameterConverter const) makeConverterTo(ParameterDefinition const & other) const = 0;

    /**
     *  @brief Interpret the given parameter vector as an ellipse, uses the given center if none is included
     *         in the parameter vector itself.
     */
    virtual afw::geom::ellipses::Ellipse makeEllipse(
        ndarray::Array<double const,1,1> const & input,
        afw::geom::Point2D const & center=afw::geom::Point2D()
    ) const = 0;

    /// Find the registered definition with the given name.
    static ParameterDefinition const & lookup(std::string const & name_);

    virtual ~ParameterDefinition() {}

protected:

    /// Construct and register a ParameterDefinition singleton instance
    explicit ParameterDefinition(std::string const & name_, int const size_);

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_parameters_h_INCLUDED
