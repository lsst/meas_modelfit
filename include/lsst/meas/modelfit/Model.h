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

#ifndef LSST_MEAS_MODELFIT_Model_h_INCLUDED
#define LSST_MEAS_MODELFIT_Model_h_INCLUDED

#include <vector>

#include "lsst/base.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/modelfit/common.h"

namespace lsst { namespace meas { namespace modelfit {

class Model;
class Prior;

struct LocalUnitTransform;

typedef std::vector<std::shared_ptr<Model>> ModelVector;

/**
 *  @brief Abstract base class and concrete factories that define multi-shapelet galaxy models
 *
 *  A Model provides a mapping from its parameters to ellipses and to a realization using
 *  shapelet objects.  A Model does not "hold" its parameters; parameters are always stored in
 *  separate arrays.
 *
 *  Model parameters are split into three categories: nonlinear, amplitudes, and fixed.  These
 *  are described more fully in @ref modelfitParameters.
 *
 *  A few private concrete subclasses of Model have been provided that will meet most needs;
 *  instances can be constructed via the make() and makeGaussian()
 */
class Model {
public:

    enum CenterEnum {
        FIXED_CENTER  = 0x0,
        SINGLE_CENTER = 0x1,
        MULTI_CENTER  = 0x2
    };

    typedef std::vector<std::string> NameVector;
    typedef std::vector<std::shared_ptr<shapelet::MultiShapeletBasis>> BasisVector;
    typedef std::vector<afw::geom::ellipses::Ellipse> EllipseVector;
    typedef std::vector<afw::geom::ellipses::Ellipse>::iterator EllipseIterator;
    typedef std::vector<afw::geom::ellipses::Ellipse>::const_iterator EllipseConstIterator;

    /**
     *  Construct a concrete Model instance with multiple ellipses and multishapelet bases
     *
     *  This can be used to construct a multi-component model (for instance, a bulge-disk decomposition),
     *  using MultiShapeletBasis objects such as those loaded by the lsst.shapelet.tractor module.
     *
     *  @param[in]  basisVector        A vector of MultiShapeletBasis objects, one for each component.
     *                                 Each component will have a separate set of ellipse parameters.
     *  @param[in]  prefixes           A vector of parameter name prefixes, one for each basis.
     *                                 These will be prepended to the names described in the documentation
     *                                 for the other overload of make().
     *  @param[in]  center             An enum specifying whether the different components should have a
     *                                 fixed center (FIXED_CENTER), the same center (SINGLE_CENTER), or
     *                                 independent centers (MULTI_CENTER).
     */
    static std::shared_ptr<Model> make(BasisVector basisVector, NameVector const & prefixes, CenterEnum center);

    /**
     *  Construct a concrete Model instance with a single ellipse and multishapelet basis
     *
     *  This can be used to construct a single-component model (e.g. a single fixed-index Sersic profile),
     *  or a linear combination model with only one ellipse.
     *
     *  @param[in]  basis              A MultiShapeletBasis object, of the sort provided by the
     *                                 lsst.shapelet.tractor module.
     *  @param[in]  center             An enum specifying whether the model should have a fixed center
     *                                 (FIXED_CENTER) or parametrized center (SINGLE_CENTER or MULTI_CENTER).
     *
     *  The names of the nonlinear and fixed parameters will be ["eta1", "eta2", "logR", "x", "y"],
     *  with "x" and "y" in the fixed parameters if center==FIXED_CENTER.  The amplitudes will be labeled
     *  "alphaN", where "N" is an integer starting from 0.
     *
     *  As implied by the parameter names, the ellipse is parametrized using
     *  afw::geom::ellipses::SeparableConformalShearLogTraceRadius.  For the basis objects provided
     *  by lsst.shapelet.tractor, that generally means that logR=0 corresponds to the half-light radius.
     */
    static std::shared_ptr<Model> make(std::shared_ptr<shapelet::MultiShapeletBasis> basis, CenterEnum center);

    /**
     *  Construct a concrete Model instance that represents a single elliptical Gaussian function.
     *
     *  The Models returned by this method use the same ellipse parametrization and naming schemes as
     *  those returned by make().
     *
     *  @param[in]  center             An enum specifying whether the model should have a fixed center
     *                                 (FIXED_CENTER) or parametrized center (SINGLE_CENTER or MULTI_CENTER).
     *  @param[in]  radius             The radius at which logR=0, in units of the Gaussian sigma parameter
     *                                 (i.e. radius=1 corresponds to a model in which the radius parameter
     *                                 is ln(sigma)).
     */
    static std::shared_ptr<Model> makeGaussian(CenterEnum center, double radius=1.0);

    /// Return the number of free nonlinear parameters
    int getNonlinearDim() const { return _nonlinearNames.size(); }

    /// Return the number of linear parameters
    int getAmplitudeDim() const { return _amplitudeNames.size(); }

    /// Return the number of fixed nonlinear parameters
    int getFixedDim() const { return _fixedNames.size(); }

    /// Return the number of MultiShapeletBasis objects (equivalently, the number of ellipses)
    int getBasisCount() const { return _basisVector.size(); }

    /// Return the names of the free nonlinear parameters
    NameVector const & getNonlinearNames() const { return _nonlinearNames; }

    /// Return the names of the amplitude parameters
    NameVector const & getAmplitudeNames() const { return _amplitudeNames; }

    /// Return the names of the fixed nonlinear parameters
    NameVector const & getFixedNames() const { return _fixedNames; }

    /// Return the MultiShapeletBasis objects that comprise the Model
    BasisVector const & getBasisVector() const { return _basisVector; }

    /// Create a MultiShapeletFunction object from a set of parameter vectors.
    shapelet::MultiShapeletFunction makeShapeletFunction(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar const,1,1> const & fixed
    ) const;

    /// Given an arbitrary prior, return one compatible with this Model or throw LogicError
    virtual std::shared_ptr<Prior> adaptPrior(std::shared_ptr<Prior> prior) const = 0;

    /**
     *  Return an uninitialized vector of afw::geom::ellipses::Ellipse with the parametrization expected
     *  by readEllipses() and writeEllipses().
     */
    virtual EllipseVector makeEllipseVector() const = 0;

    /**
     *  @brief Convert a set of nonlinear+fixed parameter arrays to a vector of ellipses.
     *
     *  @param[in] nonlinearIter    Pointer to the beginning of a nonlinear parameter array.
     *  @param[in] fixedIter        Pointer to the beginning of a fixed parameter array.
     *  @param[out] ellipseIter     Iterator to the beginning of an ellipse vector, as returned
     *                              by makeEllipseVector().
     *
     *  @warning The ellipse iterator *must* point to an EllipseVector originally constructed by
     *           makeEllipseVector()
     *
     *  @warning Calling writeEllipses() followed by readEllipses() does not guarantee that that the
     *           parameters on output will be the same as those on input, as the parameters may be
     *           degenerate.  However, calling readEllipses() followed by writeEllipses() is guaranteed
     *           to round-trip the ellipses.
     */
    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const = 0;

    /**
     *  @brief Convert a vector of ellipses to a set of nonlinear+fixed parameter arrays.
     *
     *  @param[in] ellipseIter     Iterator to the beginning of an ellipse vector, as returned
     *                              by makeEllipseVector().
     *  @param[out] nonlinearIter    Pointer to the beginning of a nonlinear parameter array.
     *  @param[out] fixedIter        Pointer to the beginning of a fixed parameter array.
     *
     *  @warning The ellipse iterator *must* point to an EllipseVector originally constructed by
     *           makeEllipseVector()
     *
     *  @warning Calling writeEllipses() followed by readEllipses() does not guarantee that that the
     *           parameters on output will be the same as those on input, as the parameters may be
     *           degenerate.  However, calling readEllipses() followed by writeEllipses() is guaranteed
     *           to round-trip the ellipses.
     */
    virtual void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const = 0;

    /**
     *  @brief Convert a set of nonlinear+fixed parameter arrays to a vector of ellipses.
     *
     *  @param[in] nonlinear        nonlinear parameter array.
     *  @param[in] fixed            fixed parameter array.
     *
     *  This is a convenient method that combines the call to makeEllipseVector() with a call to
     *  the other overload of writeEllipses(), for cases when there is no need to reuse an existing
     *  ellipse vector.
     */
    EllipseVector writeEllipses(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & fixed
    ) const;

    /**
     *  @brief Convert a vector of ellipses to a set of nonlinear+fixed parameter arrays.
     *
     *  @param[in] ellipses          An ellipse vector, as returned by makeEllipseVector().
     *  @param[out] nonlinear        Output nonlinear parameter array.
     *  @param[out] fixed            Output fixed parameter array.
     *
     *  @warning The EllipseVector must have been originally constructed by
     *           makeEllipseVector()
     */
    void readEllipses(
        EllipseVector const & ellipses,
        ndarray::Array<Scalar,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & fixed
    ) const;

    /**
     *  Transform (in-place) parameter vectors from one unit system to another.
     *
     *  The default implementation transforms nonlinear and fixed parameters by converting them to
     *  ellipses, transforming the ellipses, and converting back to parameters.  The amplitudes are
     *  simply multiplied by transform.flux.  Subclasses for which this isn't appropriate should
     *  override.
     */
    virtual void transformParameters(
        LocalUnitTransform const & transform,
        ndarray::Array<Scalar,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & fixed
    ) const;

    virtual ~Model() {}

    // No copying
    Model (const Model&) = delete;
    Model& operator=(const Model&) = delete;

    // No moving
    Model (Model&&) = delete;
    Model& operator=(Model&&) = delete;

protected:

    Model(
        BasisVector basisVector,
        NameVector nonlinearNames,
        NameVector amplitudeNames,
        NameVector fixedNames
    );

private:
    NameVector _nonlinearNames;
    NameVector _amplitudeNames;
    NameVector _fixedNames;
    BasisVector _basisVector;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_Model_h_INCLUDED
