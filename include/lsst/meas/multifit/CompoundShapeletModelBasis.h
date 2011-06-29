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

class CompoundShapeletImpl;

} // namespace detail

/**
 *  @brief An ModelBasis subclass for a multi-scale shapelet expansion.
 */
class CompoundShapeletModelBasis : public ModelBasis {
public:

    typedef boost::shared_ptr<CompoundShapeletModelBasis> Ptr;

    typedef std::vector<ShapeletModelBasis::Ptr> ComponentVector;
    
    ComponentVector extractComponents() const;

    lsst::ndarray::Array<Pixel const,2,1> getForward() const;

    lsst::ndarray::Array<Pixel const,2,1> getReverse() const;

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

    /// @brief Return a matrix of integral inner products of the basis functions.
    virtual Eigen::MatrixXd computeInnerProductMatrix(
        lsst::afw::geom::ellipses::BaseCore const & ellipse
    ) const;

    /**
     *  @brief Return a matrix of inner products that can be used to orthogonalize the basis.
     */
    Eigen::MatrixXd computeInnerProductMatrix() const;

    static Ptr load(std::string const & filename);
    void save(std::string const & filename);

    virtual ~CompoundShapeletModelBasis();

protected:

    virtual void _integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const;

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
    friend class CompoundShapeletBuilder;

    explicit CompoundShapeletModelBasis(boost::shared_ptr<detail::CompoundShapeletImpl> const & impl);

    boost::shared_ptr<detail::CompoundShapeletImpl> _impl;
};

/**
 *  @brief A builder class for CompoundShapeletModelBasis.
 */
class CompoundShapeletBuilder {
public:

    typedef std::vector<ShapeletModelBasis::Ptr> ComponentVector;

    explicit CompoundShapeletBuilder(ComponentVector const & components);

    explicit CompoundShapeletBuilder(
        ComponentVector const & components,
        lsst::ndarray::Array<Pixel const,2,1> const & forward,
        lsst::ndarray::Array<Pixel const,2,1> const & reverse
    );
    
    ComponentVector extractComponents() const;

    lsst::ndarray::Array<Pixel const,2,1> getForward() const;

    lsst::ndarray::Array<Pixel const,2,1> getReverse() const;

    /**
     *  @brief Return a matrix of inner products that can be used to orthogonalize the basis.
     */
    Eigen::MatrixXd computeInnerProductMatrix() const;

    void orthogonalize();

    void normalize();

    int getSize() const;

    void slice(int start, int stop);

    void setMapping(
        lsst::ndarray::Array<Pixel const,2,1> const & forward,
        lsst::ndarray::Array<Pixel const,2,1> const & reverse
    );

    CompoundShapeletModelBasis::Ptr build() const;

    ~CompoundShapeletBuilder();

private:
    void edit();
    boost::shared_ptr<detail::CompoundShapeletImpl> _impl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_CompoundShapeletModelBasis
