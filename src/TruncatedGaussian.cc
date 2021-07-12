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

//
// See the "Truncated Gaussian Math" Doxygen page for formulae and mathematical
// explanations of the functions implemented here.
//
// Some variable names are capitalized here in violation of LSST standards
// because it makes them better match the formulae on that page (which uses
// standard mathematical capitalizations.  I think the connection with the math
// is far more important that code standards adherence in this case.
//

#include "boost/math/special_functions/erf.hpp"
#include <memory>
#include "Eigen/Eigenvalues"
#include "Eigen/LU"

#include "ndarray/eigen.h"
#include "lsst/log/Log.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/modelfit/integrals.h"
#include "lsst/meas/modelfit/TruncatedGaussian.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

static double const THRESHOLD = 1E-15;
static double const SLN_THRESHOLD = 1E-5;

} // anonymous

// -------- Main TruncatedGaussian class --------------------------------------------------------------------

class TruncatedGaussian::Impl {
public:

    explicit Impl(int n) :
        untruncatedFraction(1.0), logPeakAmplitude(1.0), logIntegral(1.0),
        mu(n), s(n), v(n,n)
        {}

    Scalar untruncatedFraction;
    Scalar logPeakAmplitude;
    Scalar logIntegral;
    Vector mu;
    Vector s;  // H = Sigma^{-1} = V S V^T
    Matrix v;
};

TruncatedGaussian TruncatedGaussian::fromSeriesParameters(
    Scalar q0, Vector const & gradient, Matrix const & hessian
) {
    static Scalar const LN_2PI = std::log(2.0*M_PI);
    static Scalar const SQRT_PI = std::sqrt(M_PI);
    LOG_LOGGER trace4Logger = LOG_GET("TRACE4.meas.modelfit.TruncatedGaussian");
    int const n = gradient.size();
    if (hessian.rows() != n || hessian.cols() != n) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Mismatch between grad size (%d) and hessian dimensions (%d, %d)")
             % n % hessian.rows() % hessian.cols()).str()
        );
    }
    std::shared_ptr<Impl> impl = std::make_shared<Impl>(n);
    if (n == 1) {
        Scalar g = gradient[0];
        Scalar H = hessian(0,0);
        Scalar mu = -g / H;
        LOGL_DEBUG(trace4Logger, "fromSeriesParameters: 1d with H=[%g], mu=[%g]", H, mu);
        impl->mu[0] = mu;
        impl->s(0,0) = H;
        impl->v.setIdentity();
        impl->logPeakAmplitude = q0 + 0.5*g*mu;
        impl->untruncatedFraction = 0.5*boost::math::erfc(-mu*std::sqrt(H/2.0));
        impl->logIntegral = impl->logPeakAmplitude + 0.5*std::log(H) - 0.5*LN_2PI
            - std::log(impl->untruncatedFraction);
    } else if (n == 2) {
        // There are some unnecessary copies here, but they help keep the notation clean;
        // someday we'll be able to use auto and make them mostly references, but for now
        // it's too tricky to figure out the Eigen return types.
        Eigen::Vector2d g = gradient.segment<2>(0);
        Eigen::Matrix2d H = hessian.block<2,2>(0,0);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(H);
        Eigen::Vector2d s = eig.eigenvalues();
        Eigen::Matrix2d v = eig.eigenvectors();
        if (v(0,1) < 0.0) { // just enforces a convention
            v(0,1) *= -1;
            v(1,1) *= -1;
        }
        LOGL_DEBUG(trace4Logger,
            "fromSeriesParameters: 2d with H=[[%16.16g, %16.16g], [%16.16g, %16.16g]], g=[%16.16g, %16.16g]",
            H(0,0), H(0,1), H(1,0), H(1,1), g[0], g[1]
        );
        LOGL_DEBUG(trace4Logger, "fromSeriesParameters: v=[[%8.8g, %8.8g], [%8.8g, %8.8g]]",
                          v(0,0), v(0,1), v(1,0), v(1,1));
        Eigen::Vector2d mu;
        bool isSingular = (s[0] < s[1] * THRESHOLD);
        if (isSingular) {
            double solutionThreshold = SLN_THRESHOLD*std::max(std::abs(g[0]), std::abs(g[1]));
            if (!(std::abs(v.col(0).dot(g)) < solutionThreshold)) {
                LOGL_DEBUG(trace4Logger,
                    "no solution: singular matrix with s=[%g, %g], mu=[%g, %g]",
                    s[0], s[1], mu[0], mu[1]
                );
                LOGL_DEBUG(trace4Logger,
                    "|v.col(0).dot(g)=%g| > %g",
                    v.col(0).dot(g), solutionThreshold
                );
                throw LSST_EXCEPT(
                    pex::exceptions::RuntimeError,
                    "Integral diverges: H*mu = -g has no solution"
                );
            }
            if (v(1,1) < 0.0) {
                LOGL_DEBUG(trace4Logger, "unconstrained: singular matrix with s=[%g, %g], mu=[%g, %g]; v(1,1)=%g",
                                  s[0], s[1], mu[0], mu[1], v(1,1));
                throw LSST_EXCEPT(
                    pex::exceptions::RuntimeError,
                    "Integral diverges: degenerate direction is not constrained"
                );
            }
            mu = -(v.col(1).dot(g)/s[1]) * v.col(1);
            s[0] = 0.0; // better than leaving it as it was, because it might have been tiny and negative
            Scalar z = v.col(1).dot(g) / std::sqrt(2.0 * s[1]);
            // we use abs() here because we know we want a positive integral, but we don't know which of
            // v(0,0) and v(0,1) is positive and which is negative
            impl->logIntegral = q0 - std::log(
                std::abs(
                    (v(0,1)/v(0,0) - v(1,1)/v(1,0))
                    * (1.0 - z*SQRT_PI*std::exp(z*z)*boost::math::erfc(z))
                    / s[1]
                )
            );
            impl->untruncatedFraction = 0.0; // untruncated integral diverges, so ratio is 0
            impl->logPeakAmplitude = q0 + 0.5*g.dot(mu);
        } else {
            mu = -v * ((v.adjoint() * g).array() / s.array()).matrix();
            LOGL_DEBUG(trace4Logger, "fromSeriesParameters: full-rank matrix with s=[%g, %g], mu=[%g, %g]",
                              s[0], s[1], mu[0], mu[1]);
            Scalar rho = -H(0,1) / std::sqrt(H(0,0) * H(1,1));
            Scalar detH = s[0] * s[1];
            Scalar sigma00 = H(1,1) / detH;
            Scalar sigma11 = H(0,0) / detH;
            impl->logPeakAmplitude = q0 + 0.5*g.dot(mu);
            impl->untruncatedFraction = detail::bvnu(
                -mu[0]/std::sqrt(sigma00), -mu[1]/std::sqrt(sigma11), rho
            );
            impl->logIntegral = impl->logPeakAmplitude + 0.5*std::log(detH) - LN_2PI
                - std::log(impl->untruncatedFraction);
        }
        LOGL_DEBUG(trace4Logger, "fromSeriesParameters: v=[[%g, %g], [%g, %g]]",
                          v(0,0), v(0,1), v(1,0), v(1,1));
        impl->mu.head<2>() = mu;
        impl->s.head<2>() = s;
        impl->v.block<2,2>(0,0) = v;
    } else {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Greater than 2 dimensions not yet supported"
        );
    }
    LOGL_DEBUG(trace4Logger, "fromSeriesParameters: logPeakAmplitude=%g, logIntegral=%g, untruncatedFraction=%g",
                      impl->logPeakAmplitude, impl->logIntegral, impl->untruncatedFraction);
    return TruncatedGaussian(impl);
}

TruncatedGaussian TruncatedGaussian::fromStandardParameters(
    Vector const & mean, Matrix const & covariance
) {
    static Scalar const LN_2PI = std::log(2.0*M_PI);
    LOG_LOGGER trace4Logger = LOG_GET("TRACE4.meas.modelfit.TruncatedGaussian");
    int const n = mean.size();
    if (covariance.rows() != n || covariance.cols() != n) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Mismatch between mean size (%d) and covariance dimensions (%d, %d)")
             % n % covariance.rows() % covariance.cols()).str()
        );
    }
    std::shared_ptr<Impl> impl = std::make_shared<Impl>(n);
    if (n == 1) {
        Scalar mu = mean[0];
        Scalar Sigma = covariance(0,0);
        LOGL_DEBUG(trace4Logger, "fromStandardParameters: 1d with Sigma=[%g], mu=[%g]", Sigma, mu);
        impl->mu[0] = mu;
        impl->s(0,0) = 1.0/Sigma;
        impl->v.setIdentity();
        impl->untruncatedFraction = 0.5*boost::math::erfc(-mu/std::sqrt(2.0*Sigma));
        impl->logPeakAmplitude = std::log(impl->untruncatedFraction) + 0.5*std::log(Sigma) + 0.5*LN_2PI;
        impl->logIntegral = 0.0;
    } else if (n == 2) {
        // There are some unnecessary copies here, but they help keep the notation clean;
        // someday we'll be able to use auto and make them mostly references, but for now
        // it's too tricky to figure out the Eigen return types.
        Eigen::Vector2d mu = mean.segment<2>(0);
        Eigen::Matrix2d Sigma = covariance.block<2,2>(0,0);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(Sigma);
        Eigen::Vector2d s = eig.eigenvalues().array().inverse().matrix();
        Eigen::Matrix2d v = eig.eigenvectors();
        bool isSingular = !(s[1] >= s[0] * THRESHOLD);
        if (isSingular) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                "TruncatedGaussian cannot be normalized"
            );
        }
        LOGL_DEBUG(trace4Logger, "fromStandardParameters: full-rank matrix with s=[%g, %g], mu=[%g, %g]",
                          s[0], s[1], mu[0], mu[1]);
        Scalar rho = Sigma(0,1) / std::sqrt(Sigma(0,0) * Sigma(1,1));
        Scalar detSigma = 1.0 / (s[0] * s[1]);
        impl->untruncatedFraction = detail::bvnu(
            -mu[0]/std::sqrt(Sigma(0,0)), -mu[1]/std::sqrt(Sigma(1,1)), rho
        );
        impl->logPeakAmplitude = std::log(impl->untruncatedFraction) + 0.5*std::log(detSigma) + LN_2PI;
        impl->logIntegral = 0.0;
        impl->mu.head<2>() = mu;
        impl->s.head<2>() = s;
        impl->v.block<2,2>(0,0) = v;
    } else {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Greater than 2 dimensions not yet supported"
        );
    }
    LOGL_DEBUG(trace4Logger, "fromStandardParameters: logPeakAmplitude=%g, logIntegral=%g, untruncatedFraction=%g",
                      impl->logPeakAmplitude, impl->logIntegral, impl->untruncatedFraction);
    return TruncatedGaussian(impl);
}

int TruncatedGaussian::getDim() const {
    return _impl->mu.size();
}

Vector TruncatedGaussian::maximize() const {
    Vector result(_impl->mu);
    int const n = _impl->mu.size();
    int k = 0;
    for (int i = 0; i < n; ++i) {
        if (result[i] < 0.0) {
            ++k;
        }
    }
    if (n == k) {
        result.setZero();
    } else if (k > 0) {
        Eigen::VectorXi indices(n);
        for (int i = 0, j1 = 0, j2 = n - k; i < n; ++i) {
            if (result[i] < 0.0) {
                indices[i] = j2;
                ++j2;
            } else {
                indices[i] = j1;
                ++j1;
            }
        }
        Eigen::PermutationMatrix<Eigen::Dynamic> p(indices);
        // beta, G, nu: permuted versions of alpha, H, mu
        Matrix pv = p * _impl->v;
        Matrix G = pv * _impl->s.asDiagonal() * pv.adjoint();
        Vector nu = p * _impl->mu;
        Vector beta = Vector::Zero(n);
        Eigen::FullPivLU<Matrix> solver(G.topLeftCorner(n - k, n - k));
        beta.head(n - k) = solver.solve(G.topRightCorner(n - k, k) * nu.tail(k)) + nu.head(n - k);
        result = p.transpose() * beta;
    }
    return result;
}

Scalar TruncatedGaussian::getUntruncatedFraction() const {
    return _impl->untruncatedFraction;
}

Scalar TruncatedGaussian::getLogPeakAmplitude() const {
    return _impl->logPeakAmplitude;
}

Scalar TruncatedGaussian::getLogIntegral() const {
    return _impl->logIntegral;
}

TruncatedGaussian::~TruncatedGaussian() {}

// -------- LogEvaluator class ------------------------------------------------------------------------------

TruncatedGaussianLogEvaluator::TruncatedGaussianLogEvaluator(TruncatedGaussian const & parent) :
    _norm(parent._impl->logPeakAmplitude), _mu(parent._impl->mu), _workspace(_mu.size()),
    _rootH(parent._impl->s.array().sqrt().matrix().asDiagonal() * parent._impl->v.adjoint())
{}

Scalar TruncatedGaussianLogEvaluator::operator()(ndarray::Array<Scalar const,1,1> const & alpha) const {
    return (*this)(ndarray::asEigenMatrix(alpha));
}

void TruncatedGaussianLogEvaluator::operator()(
    ndarray::Array<Scalar const,2,1> const & alpha,
    ndarray::Array<Scalar,1,1> const & output
) const {
    LSST_THROW_IF_NE(
        alpha.getSize<0>(), output.getSize<0>(),
        pex::exceptions::LengthError,
        "Size of inputs (%d) does not match size of outputs (%d)"
    );
    for (int i = 0, n = alpha.getSize<0>(); i < n; ++i) {
        output[i] = (*this)(ndarray::asEigenMatrix(alpha[i]));
    }
}

// -------- Evaluator class ---------------------------------------------------------------------------------

Scalar TruncatedGaussianEvaluator::operator()(ndarray::Array<Scalar const,1,1> const & alpha) const {
    return (*this)(ndarray::asEigenMatrix(alpha));
}

void TruncatedGaussianEvaluator::operator()(
    ndarray::Array<Scalar const,2,1> const & alpha,
    ndarray::Array<Scalar,1,1> const & output
) const {
    LSST_THROW_IF_NE(
        alpha.getSize<0>(), output.getSize<0>(),
        pex::exceptions::LengthError,
        "Size of inputs (%d) does not match size of outputs (%d)"
    );
    for (int i = 0, n = alpha.getSize<0>(); i < n; ++i) {
        output[i] = (*this)(ndarray::asEigenMatrix(alpha[i]));
    }
}

// -------- Sampler classes ---------------------------------------------------------------------------------

class TruncatedGaussianSampler::Impl {
public:

    virtual Scalar apply(afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha) = 0;

    virtual ~Impl() {}
};

namespace {

class SamplerImplDWR1 : public TruncatedGaussianSampler::Impl {
public:

    SamplerImplDWR1(TruncatedGaussian const & parent, Vector const & mu, Matrix const & v, Vector const & s) :
        _mu(mu[0]), _rootSigma(std::sqrt(1.0/s[0]) * v(0,0))
        {}

    virtual Scalar apply(afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha) {
        do {
            alpha[0] = _rootSigma * rng.gaussian() + _mu;
        } while (alpha[0] < 0.0);
        return 1.0;
    }

private:
    Scalar _mu;
    Scalar _rootSigma;
};

class SamplerImplDWR : public TruncatedGaussianSampler::Impl {
public:

    SamplerImplDWR(TruncatedGaussian const & parent, Vector const & mu, Matrix const & v, Vector const & s) :
        _mu(mu), _workspace(mu.size()),
        _rootSigma(v * s.array().inverse().sqrt().matrix().asDiagonal() * v.adjoint())
        {}

    virtual Scalar apply(afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha) {
        do {
            for (int j = 0; j < _workspace.size(); ++j) {
                _workspace[j] = rng.gaussian();
            }
            ndarray::asEigenMatrix(alpha) = _rootSigma * _workspace + _mu;
        } while ((ndarray::asEigenArray(alpha) < 0.0).any());
        return 1.0;
    }

private:
    Vector _mu;
    Vector _workspace;
    Matrix _rootSigma;
};

Scalar draw1d(afw::math::Random & rng, Scalar Ap) {
    return -boost::math::erfc_inv(2.0*(1.0 - rng.uniform()*Ap)) * M_SQRT2;
}

class SamplerImplAAW1 : public TruncatedGaussianSampler::Impl {
public:

    SamplerImplAAW1(
        TruncatedGaussian const & parent, Vector const & mu, Matrix const & v, Vector const & s
    ) :
        _mu(mu[0]), _rootD(std::sqrt(1.0/s[0])),
        _A(0.5*boost::math::erfc(-_mu/(M_SQRT2*_rootD)))
        {}

    virtual Scalar apply(afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha) {
        alpha[0] = draw1d(rng, _A) * _rootD + _mu;
        return 1.0;
    }

private:
    Scalar _mu;
    Scalar _rootD;
    Scalar _A;
};

// We inherit from TruncatedGaussianSampler not just because we want to evaluate the function
// repeatedly, but also because we want to reuse some of its data members (mu, workspace) for
// our own purposes, and we can only do that via inheritance rather than containment.
class SamplerImplAAW : public TruncatedGaussianSampler::Impl, private TruncatedGaussianLogEvaluator {
public:

    SamplerImplAAW(
        TruncatedGaussian const & parent, Vector const & mu, Matrix const & v, Vector const & s
    ) :
        TruncatedGaussianLogEvaluator(parent),
        _pNorm(1.0),
        _lnAf(parent.getLogIntegral()),
        _Ap(mu.size()),
        _rootD(mu.size())
        {
            static Scalar const LOG_2PI = std::log(2.0*M_PI);
            static Scalar const MAX_NEGATIVE_SIGMA = 6.0;
            LOG_LOGGER trace4Logger = LOG_GET("TRACE4.meas.modelfit.TruncatedGaussian");
            // We start with the inverse of the diagonal of H; we'd prefer the diagonal of the inverse,
            // but H may not be invertible.  The inverse of the diagonal at H represents the width in
            // each dimension at a fixed point in all the other dimensions, so it's not quite wide
            // enough to be a good importance distribution, so we'll multiply it by 2 (just a heuristic).
            _rootD = (v * s.asDiagonal() * v.adjoint()).diagonal().array().inverse().sqrt().matrix() * 2;
            LOGL_DEBUG(trace4Logger, "AAW Sampler: rootD=[%g, %g]", _rootD[0], _rootD[1]);
            for (int j = 0; j < _Ap.size(); ++j) {
                Scalar x = _mu[j]/_rootD[j];
                if (x < -MAX_NEGATIVE_SIGMA) {
                    // if mu is more than MAX_NEGATIVE SIGMA away from the feasible region, we increase
                    // _rootD so that it's exactly MAX_NEGATIVE_SIGMA away (just for numerical reasons)
                    x = -MAX_NEGATIVE_SIGMA;
                    _rootD[j] = _mu[j] / x;
                }
                _Ap[j] = 0.5*boost::math::erfc(-x / M_SQRT2);
                _pNorm += 0.5*LOG_2PI + std::log(_rootD[j]*_Ap[j]);
            }
        }

    virtual Scalar apply(afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha) {
        for (int j = 0; j < _workspace.size(); ++j) {
            // Start by drawing truncated normal deviates without scaling and shifting by rootD, mu
            // because we'd have to undo that shift and scale to evaluate the proposal.
            alpha[j] = draw1d(rng, _Ap[j]);
        }
        Scalar logProposal = 0.5*ndarray::asEigenMatrix(alpha).squaredNorm() + _pNorm;
        // Now that we've evaluated the proposal, we apply the scaling and shifting
        ndarray::asEigenArray(alpha) *= _rootD.array();
        ndarray::asEigenMatrix(alpha) += _mu;
        // Call private Evaluator base class, divide by integral (in log space)
        Scalar logActual = (*this)(ndarray::asEigenMatrix(alpha)) - _lnAf;
        return std::exp(logProposal - logActual);
    }

private:
    Scalar _pNorm; // normalization factor for full N-d importance distribution
    Scalar _lnAf; // log integral of the true N-d distribution
    Vector _Ap; // untruncated fractions for each 1-d importance distribution
    Vector _rootD; // sqrt of variances for importance distribution
};

} // anonymous

TruncatedGaussianSampler::TruncatedGaussianSampler(
    TruncatedGaussian const & parent,
    TruncatedGaussian::SampleStrategy strategy
) {
    LOG_LOGGER trace4Logger = LOG_GET("TRACE4.meas.modelfit.TruncatedGaussian");
    if (parent.getDim() == 1) {
        switch (strategy) {
        case TruncatedGaussian::DIRECT_WITH_REJECTION:
            _impl = std::make_shared<SamplerImplDWR1>(
                parent, parent._impl->mu, parent._impl->v, parent._impl->s
            );
            LOGL_DEBUG(trace4Logger, "Sampler: using DWR1");
            break;
        case TruncatedGaussian::ALIGN_AND_WEIGHT:
            _impl = std::make_shared<SamplerImplAAW1>(
                parent, parent._impl->mu, parent._impl->v, parent._impl->s
            );
            LOGL_DEBUG(trace4Logger, "Sampler: using AAW1");
            break;
        }
    } else {
        switch (strategy) {
        case TruncatedGaussian::DIRECT_WITH_REJECTION:
            _impl = std::make_shared<SamplerImplDWR>(
                parent, parent._impl->mu, parent._impl->v, parent._impl->s
            );
            LOGL_DEBUG(trace4Logger, "Sampler: using DWR");
            break;
        case TruncatedGaussian::ALIGN_AND_WEIGHT:
            _impl = std::make_shared<SamplerImplAAW>(
                parent, parent._impl->mu, parent._impl->v, parent._impl->s
            );
            LOGL_DEBUG(trace4Logger, "Sampler: using AAW");
            break;
        }
    }
    if (!_impl) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Invalid enum value for SampleStrategy"
        );
    }
}

Scalar TruncatedGaussianSampler::operator()(
    afw::math::Random & rng, ndarray::Array<Scalar,1,1> const & alpha
) const {
    return _impl->apply(rng, alpha);
}

void TruncatedGaussianSampler::operator()(
    afw::math::Random & rng,
    ndarray::Array<Scalar,2,1> const & alpha,
    ndarray::Array<Scalar,1,1> const & weights,
    bool multiplyWeights
) const {
    LSST_THROW_IF_NE(
        alpha.getSize<0>(), weights.getSize<0>(),
        pex::exceptions::LengthError,
        "First dimension of alpha array (%d) does not match size of weights array (%d)"
    );
    if (multiplyWeights) {
        for (int i = 0, n = alpha.getSize<0>(); i < n; ++i) {
            weights[i] *= _impl->apply(rng, alpha[i]);
        }
    } else {
        for (int i = 0, n = alpha.getSize<0>(); i < n; ++i) {
            weights[i] = _impl->apply(rng, alpha[i]);
        }
    }
}

TruncatedGaussianSampler::~TruncatedGaussianSampler() {} // defined in .cc so it can see Impl's dtor

}}} // namespace lsst::meas::modelfit
