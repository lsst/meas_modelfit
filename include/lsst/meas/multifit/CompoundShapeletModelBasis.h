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

#ifndef LSST_MEAS_MULTIFIT_SHAPELETS_CompoundShapeletModelBasis
#define LSST_MEAS_MULTIFIT_SHAPELETS_CompoundShapeletModelBasis

#include "lsst/meas/multifit/ShapeletModelBasis.h"
#include "lsst/ndarray/eigen.h"
#include <vector>
#include <string>

namespace lsst { namespace meas { namespace multifit {

namespace detail {

template <typename Component> class CompoundShapeletImpl;

} // namespace detail

/**
 *  @brief An ModelBasis subclass for a multi-scale shapelet expansion.
 */
class CompoundShapeletModelBasis : public ModelBasis {
    typedef detail::CompoundShapeletImpl<ShapeletModelBasis> Impl;
public:

    typedef boost::shared_ptr<CompoundShapeletModelBasis> Ptr;
    typedef std::vector<ShapeletModelBasis::Ptr> ComponentVector;

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
    
    ComponentVector extractComponents() const;

    lsst::ndarray::Array<Pixel const,2,1> getMapping() const;

    Eigen::MatrixXd computeInnerProductMatrix() const;

    static Ptr load(std::string const & filename);
    void save(std::string const & filename);

    virtual ~CompoundShapeletModelBasis();

protected:

    virtual void _evaluate(
        ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    virtual void _evaluateRadialProfile(
        ndarray::Array<Pixel,2,1> const & profile,
        ndarray::Array<Pixel const,1,1> const & radii
    ) const;

private:

    friend class CompoundShapeletBuilder;

    explicit CompoundShapeletModelBasis(
        boost::shared_ptr<Impl> const & impl,
        ndarray::Array<Pixel const,2,2> const & multipoleMatrix,
        ndarray::Array<Pixel const,2,1> const & constraintMatrix,
        ndarray::Array<Pixel const,1,1> const & constraintVector
    );

    boost::shared_ptr<Impl> _impl;
};

/**
 *  @brief A simple 1d function class for use with CompoundShapeletBuilder::approximate.
 */
class ProfileFunction {
public:
    typedef boost::shared_ptr<ProfileFunction> Ptr;

    virtual double operator()(double radius) const = 0;
    virtual ~ProfileFunction() {}

    /// @brief A truncated de Vaucouleur profile (the SDSS prescription) to use with approximate().
    static Ptr makeTruncatedDeVaucouleur();

    /// @brief A truncated exponential profile (the SDSS prescription) to use with approximate().
    static Ptr makeTruncatedExponential();
};

/**
 *  @brief A builder class for CompoundShapeletModelBasis.
 */
class CompoundShapeletBuilder {
    typedef detail::CompoundShapeletImpl<ShapeletModelBasis> Impl;
public:

    typedef std::vector<ShapeletModelBasis::Ptr> ComponentVector;

    explicit CompoundShapeletBuilder(
        ComponentVector const & components,
        lsst::afw::math::shapelets::BasisTypeEnum basisType = lsst::afw::math::shapelets::HERMITE,
        bool radialOnly = false
    );

    explicit CompoundShapeletBuilder(
        ComponentVector const & components,
        lsst::ndarray::Array<Pixel const,2,1> const & mapping
    );

    /**
     *  @brief Create a builder with one basis function that approximates the given profile.
     */
    static CompoundShapeletBuilder approximate(
        ProfileFunction const & profile,
        ComponentVector const & components,
        double sersicRadius,
        double maxRadius,
        lsst::ndarray::Array<Pixel const,1,1> const & matchRadii
    );

    int getSize() const;

    /// @brief Normalize the basis so the flux of the nth basis function is one.
    void normalizeFlux(int n);

    void orthogonalize();

    void slice(int start, int stop);
    
    ComponentVector extractComponents() const;

    lsst::ndarray::Array<Pixel const,2,1> getMapping() const;

    void setMapping(lsst::ndarray::Array<Pixel const,2,1> const & mapping);

    Eigen::MatrixXd computeInnerProductMatrix() const;

    /// @brief Evaluate the basis functions on the given footprint.
    void evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    /// @brief Fill a matrix that evaluates the radial profile of a basis expansion.
    void evaluateRadialProfile(
        lsst::ndarray::Array<Pixel,2,1> const & profile,
        lsst::ndarray::Array<Pixel const,1,1> const & radii
    ) const;

    void evaluateIntegration(lsst::ndarray::Array<Pixel,1,1> const & vector) const;

    void evaluateMultipoleMatrix(lsst::ndarray::Array<Pixel,2,2> const & matrix) const;

    /// @brief Return the number of inequality constraints.
    int getConstraintCount() const { return _constraintMatrix.getSize<0>(); }

    lsst::ndarray::Array<Pixel const,2,1> getConstraintMatrix() const { return _constraintMatrix; }
    lsst::ndarray::Array<Pixel const,1,1> getConstraintVector() const { return _constraintVector; }

    void setConstraint(
        lsst::ndarray::Array<Pixel const,2,1> const & matrix,
        lsst::ndarray::Array<Pixel const,1,1> const & vector
    );

    CompoundShapeletModelBasis::Ptr build() const;

    ~CompoundShapeletBuilder();

private:

    explicit CompoundShapeletBuilder(boost::shared_ptr<Impl> const & impl) : _impl(impl) {}

    boost::shared_ptr<Impl> _impl;
    ndarray::Array<Pixel,2,1> _constraintMatrix;
    ndarray::Array<Pixel,1,1> _constraintVector;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_CompoundShapeletModelBasis
