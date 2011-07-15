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

class CompoundShapeletBase {
public:

    typedef std::vector<ShapeletModelBasis::Ptr> ComponentVector;
    
    ComponentVector extractComponents() const;

    /// @brief Return the mapping from Shapelet components (rows) to compound basis functions (columns).
    lsst::ndarray::Array<Pixel const,2,1> getMapping() const { return _mapping; }

    /**
     *  @brief Return a matrix of inner products that can be used to orthogonalize the basis.
     */
    Eigen::MatrixXd computeInnerProductMatrix() const;

    void integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const;

protected:

    typedef ndarray::EigenView<Pixel,2,1> Matrix;

    struct Element {
        ShapeletModelBasis::Ptr component;
        Matrix mapping;

        Element(
            ShapeletModelBasis::Ptr const & component_, 
            ndarray::Array<Pixel,2,1> const & fullMapping,
            int offset
        );

        Element & operator=(Element const & other);
    };

    typedef std::vector<Element> ElementVector;

    CompoundShapeletBase(
        ComponentVector const & components,
        ndarray::Array<Pixel,2,1> const & mapping
    );
    
    explicit CompoundShapeletBase(ComponentVector const & components);
    void _fillElements(ComponentVector const & components);

    void _resetElements();

    static int _computeSize(ComponentVector const & components);

    static ndarray::Array<Pixel,2,2> _makeIdentity(int size);

    ElementVector _elements;
    ndarray::Array<Pixel,2,1> _mapping;
};

} // namespace detail

class CompoundShapeletBuilder;

/**
 *  @brief An ModelBasis subclass for a multi-scale shapelet expansion.
 */
class CompoundShapeletModelBasis : public ModelBasis, public detail::CompoundShapeletBase {
public:

    typedef boost::shared_ptr<CompoundShapeletModelBasis> Ptr;

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

    static Ptr load(std::string const & filename);
    void save(std::string const & filename);

    void integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
        detail::CompoundShapeletBase::integrate(vector);
    }

    virtual ~CompoundShapeletModelBasis() {}
protected:

    virtual void _integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
        detail::CompoundShapeletBase::integrate(vector);
    }

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
    FRIEND_MAKE_SHARED_1(CompoundShapeletModelBasis, lsst::meas::multifit::CompoundShapeletBuilder);

    explicit CompoundShapeletModelBasis(CompoundShapeletBuilder const &);
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
class CompoundShapeletBuilder : public detail::CompoundShapeletBase {
public:

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

    /// @brief Normalize the basis so the flux of the nth basis function is one.
    void normalizeFlux(int n);

    void orthogonalize();

    void slice(int start, int stop);

    int getSize() const { return _mapping.getSize<1>(); }

    void setMapping(
        lsst::ndarray::Array<Pixel const,2,1> const & mapping
    );

    void setConstraint(
        lsst::ndarray::Array<Pixel const,2,1> const & matrix,
        lsst::ndarray::Array<Pixel const,1,1> const & vector
    );

    CompoundShapeletModelBasis::Ptr build() const;

    lsst::ndarray::Array<Pixel,2,2> dataImage;
    lsst::ndarray::Array<Pixel,3,3> modelImages;

private:

    friend class CompoundShapeletModelBasis;

    ndarray::Array<Pixel,2,1> _constraintMatrix;
    ndarray::Array<Pixel,1,1> _constraintVector;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_CompoundShapeletModelBasis
