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

#ifndef LSST_MEAS_MULTIFIT_ModelBasis
#define LSST_MEAS_MULTIFIT_ModelBasis

#include "lsst/meas/multifit/constants.h"
#include <boost/noncopyable.hpp>

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief An abstract base class for parametrized sets of basis functions.
 *
 *  ModelBasis subclasses should be immutable.
 */
class ModelBasis : private boost::noncopyable {
public:

    typedef boost::shared_ptr<ModelBasis> Ptr;

    /**
     *  @brief Convolve the basis with the given local PSF, returning a new basis with the same
     *         parametrization.
     */
    virtual ModelBasis::Ptr convolve(CONST_PTR(LocalPsf) const & psf) const;

    /// @brief Number of basis functions.
    int getSize() const { return _size; };

    /// @brief Evaluate the basis functions on the given footprint.
    void evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    virtual ~ModelBasis() {}

protected:

    explicit ModelBasis(int size) : _size(size) {}

    ModelBasis(ModelBasis const & other) : _size(other._size) {}

    virtual void _evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const = 0;

private:

    int const _size;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_ModelBasis
