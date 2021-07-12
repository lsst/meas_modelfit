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

#ifndef LSST_MEAS_MODELFIT_Likelihood_h_INCLUDED
#define LSST_MEAS_MODELFIT_Likelihood_h_INCLUDED

#include "ndarray_fwd.h"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/modelfit/common.h"
#include "lsst/meas/modelfit/Model.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  @brief Base class for optimizer/sampler likelihood functions that compute likelihood at a point.
 *
 *  Likelihood abstracts the problem of computing the likelihood over different kinds of data.  It is
 *  responsible for creating a "model matrix" that maps amplitudes to data values, and maintaining a
 *  vector of scaled, weighted data values that corresponds to it.  Its components can be represented
 *  best in the mathematical formula for a -log likelihood assuming Gaussian data and a model with
 *  both nonlinear parameters @f$\theta@f$ and linear ("amplitude") parameters @f$\alpha@f$:
 *  @f[
 *    L(\alpha,\theta) = \frac{1}{2}\left(y - A(\theta)\alpha\right)^T\,
 *                                  \Sigma^{-1}\,\left(y - A(\theta)\alpha\right)
 *  @f]
 *  where @f$y@f$ is the data vector, @f$\Sigma@f$ is the data covariance matrix (assumed to be diagonal),
 *  and @f$A(\theta)@f$ is the "true" model matrix (parametrized on the nonlinear parameters).
 *
 *  When fitting or sampling from the likelihood, however, we don't want to use these quantities directly,
 *  and they aren't what the Likelihood class provides.  Instead, we reparametrize with:
 *  @f[
 *    w_i \equiv \Sigma_{i,i}^{-1/2}
 *  @f]
 *  @f[
 *    z_i = w_i y_i
 *  @f]
 *  @f[
 *    B_{i,j} = w_i A_{i,j}
 *  @f]
 *  resulting in the equivalent formula:
 *  @f[
 *    L(\alpha,\theta) = \frac{1}{2}\left(z-B(\theta)\alpha\right)^T\,\left(z-B(\theta)\alpha\right)
 *  @f]
 *  The @f$w_i@f$ are the weights, which are applied to both the data vector and the model matrix
 *  to account for the noise in the data.  In some cases, we may choose to use a constant weight
 *  rather than per-pixel weights, but will will still use a vector to represent it.
 */
class Likelihood
{
public:

    /// Return the number of data points
    int getDataDim() const { return _data.getSize<0>(); }

    /// Return the number of linear parameters (columns of the model matrix)
    int getAmplitudeDim() const { return _model->getAmplitudeDim(); }

    /// Return the number of nonlinear parameters (which parameterize the model matrix)
    int getNonlinearDim() const { return _model->getNonlinearDim(); }

    /// Return the number of fixed nonlinear parameters (set on Likelihood construction)
    int getFixedDim() const { return _model->getFixedDim(); }

    /// Return the vector of fixed nonlinear parameters
    ndarray::Array<Scalar const,1,1> getFixed() const { return _fixed; }

    /// Return the vector of weighted, scaled data points @f$z@f$
    ndarray::Array<Pixel const,1,1> getData() const { return _data; }

    /// Return the vector of unweighted data points @f$y@f$
    ndarray::Array<Pixel const,1,1> getUnweightedData() const { return _unweightedData; }

    /**
     *  Return the vector of weights @f$w@f$ applied to data points and model matrix rows
     *
     *  Will be an empty array if no weights are applied.
     */
    ndarray::Array<Pixel const,1,1> getWeights() const { return _weights; }

    /// Return the vector of per-data-point variances.
    ndarray::Array<Pixel const,1,1> getVariance() const { return _variance; }

    /// Return an object that defines the model and its parameters
    std::shared_ptr<Model> getModel() const { return _model; }

    /**
     *  @brief Evaluate the model for the given vector of nonlinear parameters.
     *
     *  @param[out] modelMatrix  The dataDim x amplitudeDim matrix @f$B@f$ that expresses the model
     *                           projected in such a way that it can be compared to the data
     *                           when multiplied by an amplitude vector @f$\alpha@f$.
     *                           It should be weighted if the data vector is.  The caller
     *                           is responsible for guaranteeing that the shape of the matrix
     *                           correct, but implementations should not assume anything about
     *                           the initial values of the matrix elements.
     *  @param[in] nonlinear     Vector of nonlinear parameters at which to evaluate the model.
     *  @param[in] doApplyWeights   If False, do not apply the weights to the modelMatrix.
     */
    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        bool doApplyWeights=true
    ) const = 0;

    virtual ~Likelihood() {}

    // No copying
    Likelihood ( const Likelihood & ) = delete;
    Likelihood & operator= ( const Likelihood & ) = delete;

    // No moving
    Likelihood ( Likelihood && ) = delete;
    Likelihood & operator= ( Likelihood && ) = delete;

protected:

    Likelihood(std::shared_ptr<Model> model, ndarray::Array<Scalar const,1,1> const & fixed) :
        _model(model), _fixed(fixed) {
        LSST_THROW_IF_NE(
            fixed.getSize<0>(), static_cast<std::size_t>(model->getFixedDim()),
            pex::exceptions::LengthError,
            "Fixed parameter vector size (%d) does not match Model fixed parameter dimensionality (%d)"
        );
    }

    std::shared_ptr<Model> _model;
    ndarray::Array<Scalar const,1,1> _fixed;
    ndarray::Array<Pixel,1,1> _data;
    ndarray::Array<Pixel,1,1> _unweightedData;
    ndarray::Array<Pixel,1,1> _variance;
    ndarray::Array<Pixel,1,1> _weights;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_Likelihood_h_INCLUDED
