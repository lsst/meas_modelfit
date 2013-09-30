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

#ifndef LSST_MEAS_MULTIFIT_Likelihood_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Likelihood_h_INCLUDED

#include "ndarray_fwd.h"

#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/models.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for optimizer/sampler likelihood functions that compute likelihood at a point.
 *
 *  Likelihood abstracts the problem of computing the likelihood over different kinds of data; we
 *  expect to have two major subclasses: one for single-epoch modeling, and one for multi-epoch
 *  modeling.
 */
class Likelihood
#ifndef SWIG
 : private boost::noncopyable
#endif
{
public:

    /// Return the number of data points
    int getDataDim() const { return _data.getSize<0>(); }

    /// Return the number of linear parameters (columns of the model matrix)
    int getAmplitudeDim() const { return _model->getAmplitudeDim(); }

    /// Return the number of nonlinear parameters (which parameterize the model matrix)
    int getParameterDim() const { return _model->getParameterDim(); }

    /// Return the vector of data points (possibly weighted)
    ndarray::Array<Pixel const,1,1> getData() const { return _data; }

    /// Return an object that defines the model and its parameters
    PTR(Model) getModel() const { return _model; }

    /**
     *  @brief Evaluate the model for the given vector of nonlinear parameters.
     *
     *  @param[out] modelMatrix  The dataDim x amplitudeDim matrix that expresses the model
     *                           projected in such a way that it can be compared to the data
     *                           when multiplied by an amplitude vector.
     *                           It should be weighted if the data vector is.  The caller
     *                           is responsible for guaranteeing that the shape of the matrix
     *                           correct, but implementations should not assume anything about
     *                           the initial values of the matrix elements.
     *  @param[in] parameters    Vector of nonlinear parameters at which to evaluate the model.
     */
    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & parameters
    ) const = 0;

    virtual ~Likelihood() {}

protected:

    explicit Likelihood(PTR(Model) model) : _model(model) {}

    PTR(Model) _model;
    ndarray::Array<Pixel,1,1> _data;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Likelihood_h_INCLUDED
