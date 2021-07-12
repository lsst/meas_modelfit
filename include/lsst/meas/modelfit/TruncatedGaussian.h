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

#ifndef LSST_MEAS_MODELFIT_TruncatedGaussian_h_INCLUDED
#define LSST_MEAS_MODELFIT_TruncatedGaussian_h_INCLUDED

#include "Eigen/Core"
#include "ndarray.h"

#include "lsst/base.h"
#include "lsst/afw/math/Random.h"
#include "lsst/meas/modelfit/common.h"

// TODO: we should really integrate this with Mixture somehow

namespace lsst { namespace meas { namespace modelfit {

class TruncatedGaussianSampler;
class TruncatedGaussianEvaluator;
class TruncatedGaussianLogEvaluator;

/**
 *  @brief Represents a multidimensional Gaussian function truncated at zero
 *
 *  This is typically used to represent the posterior probability of amplitude
 *  parameters, given a flat prior; we require that the amplitudes each be positive
 *  but otherwise do not modify the Gaussian likelihood.
 *
 *  Currently only 1 and 2 dimensions are supported, and all dimensions must be truncated.
 *  Computing integrals is the only operation for which > 2-d is not implemented,
 *  but the integrals must be computed upon construction, so we can't support any other
 *  operations for > 2-d either.
 *
 *  Many operations on TruncatedGaussians are defined in -log space, as underflow/overflow
 *  problems will often occur in the non-log forms.
 *
 *  See @ref modelfitTruncatedGaussianMath for implementation notes
 */
class TruncatedGaussian {
public:

    /// Enum that describes different ways of sampling from a multidimensional TruncatedGaussian
    enum SampleStrategy {

        DIRECT_WITH_REJECTION, /**< Draw from the untruncated Gaussian, and discard negative draws.
                                *   This approach is fast for minimally-truncated Gaussians, and produces
                                *   ideally distributed samples with unit weights (because we discard with
                                *   replacement).  This approach can be extremely inefficient when the
                                *   untruncated fraction is small, and cannot be used when the Hessian matrix
                                *   is singular (as this sets the untruncated fraction to zero).
                                */

        ALIGN_AND_WEIGHT       /**< Create a similar Gaussian with no x-y covariance, and importance sample
                                *   by drawing from the independent 1-d truncated Gaussians.
                                *   This approach is faster for Gaussians that are heavily truncated, but
                                *   it requires the use of per-sample weights and may produce "noisier"
                                *   samples (in that regions with low weight may have many samples, while
                                *   regions with high weight have few).  It may also be slower for
                                *   minimally truncated Gaussians because each 1-d variate requires
                                *   an inverse erfc evaluation.
                                *   In 1-d, this produces unit weights and ideally distributed samples,
                                *   but still uses an inverse erfc approach rather than rejection to
                                *   create samples.
                                */
    };

    typedef TruncatedGaussianSampler Sampler;
    typedef TruncatedGaussianLogEvaluator LogEvaluator;
    typedef TruncatedGaussianEvaluator Evaluator;

    /**
     *  @brief Create from the first and second logarithmic derivatives of the Gaussian
     *
     *  This causes the TruncatedGaussian to model the function:
     *  @f[
     *    f(\alpha) = 1_{\alpha \ge 0}(\alpha) \; e^{-q(0) -g^T \alpha -\frac{1}{2}\alpha^T H \alpha}
     *              = 1_{\alpha \ge 0}(\alpha) \; e^{-q(\alpha)}
     *  @f]
     *  Where @f$1_{\alpha \ge 0}(\alpha)@f$ is the indicator function
     *  @f[
     *    1_{\alpha \ge 0}(\alpha) \equiv \begin{cases}
     *       1 & \text{if $\alpha \ge 0$}\\
     *       0 & \text{if $\alpha \lt 0$}
     *    \end{cases}
     *  @f]
     *  Note that this implies
     *  @f[
     *    \left.\frac{\partial q}{\partial \alpha} \right|_{\alpha=0} = g
     *  @f]
     *  @f[
     *    \left.\frac{\partial^2 q}{\partial \alpha^2} \right|_{\alpha=0} = H
     *  @f]
     *  with (for @f$\alpha \ge 0@f$)
     *  @f[
     *    q(\alpha) \equiv -\ln f(\alpha)
     *  @f]
     *
     *  Note that unlike that constructed by fromStandardParameters(), in this case @f$f(\alpha)@f$ is
     *  NOT normalized to integrate to one.
     *
     *  @param[in] q0           Zeroth-order term in the expansion
     *  @param[in] gradient     Vector of first derivatives @f$g@f$
     *  @param[in] hessian      Matrix of second derivatives @f$H@f$ (symmetric)
     */
    static TruncatedGaussian fromSeriesParameters(Scalar q0, Vector const & gradient, Matrix const & hessian);

    /**
     *  @brief Create from the "standard" mean and covariance parameters of the normal distribution
     *
     *  This causes the TruncatedGaussian to model the function:
     *  @f[
     *    f(\alpha) = \frac{1_{\alpha \ge 0}(\alpha)}{k\left|2\pi\Sigma\right|^{1/2}}
     *      e^{-\frac{1}{2}(\alpha-\mu)^T \Sigma^{-1} (\alpha-\mu)}
     *    = \frac{1_{\alpha \ge 0}(\alpha)}{k}\;\Phi(\alpha-\mu,\Sigma)
     *  @f]
     *  Where @f$1_{\alpha \ge 0}(\alpha)@f$ is the indicator function
     *  @f[
     *    1_{\alpha \ge 0}(\alpha) \equiv \begin{cases}
     *       1 & \text{if $\alpha \ge 0$}\\
     *       0 & \text{if $\alpha \lt 0$}
     *    \end{cases}
     *  @f]
     *  and @f$k@f$ is the normalization constant that @f$f(\alpha)@f$ integrates to unity:
     *  @f[
     *    k \equiv \int_0^{\infty}\!\!\Phi(\alpha-\mu,\Sigma) d\alpha
     *  @f]
     *
     *  Note that unlike that constructed by fromSeriesParameters(), in this case @f$f(\alpha@f$ is
     *  normalized to integrate to one.
     *
     *  @param[in] mean       Mean vector @f$\mu@f$
     *  @param[in] covariance Covariance matrix @f$\Sigma@f$ (symmetric)
     */
    static TruncatedGaussian fromStandardParameters(Vector const & mean, Matrix const & covariance);

    /**
     *  @brief Create a Sampler object that uses the given strategy
     *
     *  @warning do not call this method with DIRECT_WITH_REJECTION unless you are certain
     *  the Gaussian is far from degenerate, as an endless loop may result in this case.
     *  It is generally safer to use the overload that selects the strategy automatically.
     */
    Sampler sample(SampleStrategy strategy) const;

    /**
     *  @brief Create a Sampler object that determines the strategy to use automatically
     *
     *  If the efficiency (number of accepted samples divided by number of total samples)
     *  is would (on average) be greater than the given value, DIRECT_WITH_REJECTION is used.
     *  If not, ALIGN_AND_WEIGHT is used.  Note that the sampling effeciency for
     *  DIRECT_WITH_REJECTION is exactly the same as the "untruncated fraction".
     */
    Sampler sample(Scalar minRejectionEfficiency=0.1) const;

    /// @brief Create a LogEvaluator object that can be used to efficiently evaluate the -log of the function
    LogEvaluator evaluateLog() const;

    /// @brief Create an Evaluator object that can be used to efficiently evaluate the function
    Evaluator evaluate() const;

    /// @brief Return the dimensionality of the function
    int getDim() const;

    /**
     *  @brief Return the location of the maximum of the truncated Gaussian.
     *
     *  This is simply the untruncated location parameter mu if all of its elements are positive;
     *  otherwise, one or more elements will be zero (and the rest may not be same as the elements
     *  of mu).
     */
    Vector maximize() const;

    /**
     *  @brief Return the fraction of the Gaussian integral that was truncated by the bounds
     *
     *  In series form (see fromSeriesParameters), this is
     *  @f[
     *    \frac{\int_0^{\infty} e^{-q(\alpha)} d\alpha}{\int_{-infty}^{\infty} e^{-q(\alpha)}  d\alpha}
     *  @f]
     *  while in standard parameter form it is simply the normalization @f$k@f$ described
     *  in fromStandardParameters().
     */
    Scalar getUntruncatedFraction() const;

    /**
     *  @brief Return the -log of the peak amplitude of the untruncated function
     *
     *  In series form, this is equal to @f$q(\mu)@f$ where @f$\mu=-H^{-1}g@f$.
     *  In standard parameter form, this is @f$\frac{1}{2}\ln\left|2\pi\Sigma\right| - \ln k@f$.
     */
    Scalar getLogPeakAmplitude() const;

    /**
     *  @brief Return the -log of the integral of the truncated function
     *
     *  More precisely, this is simply
     *  @f[
     *    -\ln\int_{-\infty}^{\infty} f(\alpha) d\alpha
     *  @f]
     *  In series form, this is equivalent to
     *  @f[
     *    -\ln\int_0^{\infty} e^{-q(\alpha)} d\alpha
     *  @f]
     *  In standard parameter form it is 0 by construction.
     */
    Scalar getLogIntegral() const;

    ~TruncatedGaussian();

private:

    friend class TruncatedGaussianSampler;
    friend class TruncatedGaussianLogEvaluator;

    class Impl;

    explicit TruncatedGaussian(std::shared_ptr<Impl> impl) : _impl(impl) {}

    std::shared_ptr<Impl> _impl;
};

/**
 *  @brief Helper class for evaluating the -log of a TruncatedGaussian
 */
class TruncatedGaussianLogEvaluator {
public:

    explicit TruncatedGaussianLogEvaluator(TruncatedGaussian const & parent);

    template <typename Derived>
    Scalar operator()(Eigen::MatrixBase<Derived> const & alpha) const {
        if ((alpha.array() < 0.0).any()) return std::numeric_limits<Scalar>::infinity();
        _workspace = alpha - _mu;
        return 0.5*(_rootH*_workspace).squaredNorm() + _norm;
    }

    Scalar operator()(ndarray::Array<Scalar const,1,1> const & alpha) const;

    void operator()(
        ndarray::Array<Scalar const,2,1> const & alpha,
        ndarray::Array<Scalar,1,1> const & output
    ) const;

protected:
    Scalar _norm;
    Vector _mu;
    mutable Vector _workspace;
    Matrix _rootH;
};

/**
 *  @brief Helper class for evaluating the -log of a TruncatedGaussian
 */
class TruncatedGaussianEvaluator {
public:

    explicit TruncatedGaussianEvaluator(TruncatedGaussian const & parent) :
        _internal(parent)
    {}

    template <typename Derived>
    Scalar operator()(Eigen::MatrixBase<Derived> const & alpha) const {
        return std::exp(-_internal(alpha));
    }

    Scalar operator()(ndarray::Array<Scalar const,1,1> const & alpha) const;

    void operator()(
        ndarray::Array<Scalar const,2,1> const & alpha,
        ndarray::Array<Scalar,1,1> const & output
    ) const;

private:
    TruncatedGaussianLogEvaluator _internal;
};

/**
 *  @brief Helper class for drawing samples from a TruncatedGaussian
 */
class TruncatedGaussianSampler {
public:

    explicit TruncatedGaussianSampler(
        TruncatedGaussian const & parent,
        TruncatedGaussian::SampleStrategy strategy
    );

    /**
     *  @brief Draw a single sample from a TruncatedGaussian
     *
     *  @param[in]  rng      Random number generator
     *  @param[out] alpha    Output sample vector to fill
     *
     *  @return the weight of the sample (always betweeen 0 and 1)
     */
    Scalar operator()(afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha) const;

    /**
     *  @brief Draw multiple samples from a TruncatedGaussian
     *
     *  @param[in]  rng      Random number generator
     *  @param[out] alpha    Output sample vector to fill; first dimension sets the number of samples
     *  @param[out] weights  Output weight vector to fill
     *  @param[in]  multiplyWeights  If true, multiply the weights vector by the weights rather than
     *                               fill it.
     */
    void operator()(
        afw::math::Random & rng,
        ndarray::Array<Scalar,2,1> const & alpha,
        ndarray::Array<Scalar,1,1> const & weights,
        bool multiplyWeights=false
    ) const;

    ~TruncatedGaussianSampler(); // defined in .cc so it can see Impl's dtor

    class Impl; // public so we can inherit from it in the .cc file

private:
    std::shared_ptr<Impl> _impl;
};

inline TruncatedGaussian::Sampler TruncatedGaussian::sample(SampleStrategy strategy) const {
    return Sampler(*this, strategy);
}

inline TruncatedGaussian::Sampler TruncatedGaussian::sample(Scalar minRejectionEfficiency) const {
    return Sampler(
        *this,
        (getUntruncatedFraction() < minRejectionEfficiency) ? ALIGN_AND_WEIGHT : DIRECT_WITH_REJECTION
    );
}

inline TruncatedGaussian::LogEvaluator TruncatedGaussian::evaluateLog() const { return LogEvaluator(*this); }

inline TruncatedGaussian::Evaluator TruncatedGaussian::evaluate() const { return Evaluator(*this); }

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_TruncatedGaussian_h_INCLUDED
