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

#ifndef LSST_MEAS_MULTIFIT_psf_h_INCLUDED
#define LSST_MEAS_MULTIFIT_psf_h_INCLUDED

#include "boost/scoped_ptr.hpp"

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/Likelihood.h"

namespace lsst { namespace meas { namespace multifit {

PTR(Model) makeMultiShapeletPsfModel(std::vector<int> const & orders);

class MultiShapeletPsfLikelihood : public Likelihood {
public:

    MultiShapeletPsfLikelihood(
        ndarray::Array<Pixel const,2,2> const & image,
        afw::geom::Point2I const & xy0,
        PTR(Model) model,
        Scalar sigma,
        ndarray::Array<Scalar const,1,1> const & fixed
    );

    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        bool doApplyWeights=true
    ) const;

    virtual ~MultiShapeletPsfLikelihood();

private:
    class Impl;
    boost::scoped_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_psf_h_INCLUDED
