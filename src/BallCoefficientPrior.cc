#include "lsst/meas/multifit/BallCoefficientPrior.h"
#include "lsst/ndarray/eigen.h"
#include "Eigen/SVD"
#include "boost/math/special_functions/gamma.hpp"

namespace lsst { namespace meas { namespace multifit {

namespace {

template <typename RowT, typename VectorT>
void drawDiagonalGaussian(
    Random & engine,
    Eigen::MatrixBase<RowT> const & x,
    Eigen::MatrixBase<VectorT> const & s,
    Eigen::MatrixBase<VectorT> const & nu
) {
    Eigen::MatrixBase<RowT> & y = const_cast<Eigen::MatrixBase<RowT> &>(x);
    for (int j = 0; j < x.size(); ++j) {
        y[j] = engine.gaussian();
    }
    y.cwise() /= s;
    y += nu;
}

template <typename RowT>
void drawUniformOnSphere(
    Random & engine,
    Eigen::MatrixBase<RowT> const & x, 
    double radius
) {
    Eigen::MatrixBase<RowT> & y = const_cast<Eigen::MatrixBase<RowT> &>(x);
    for (int j = 0; j < x.size(); ++j) {
        y[j] = engine.gaussian();
    }
    y *= std::pow(engine.uniform(), 1.0 / x.size()) / x.norm() * radius;
}

} // anonymous

struct BallCoefficientPrior::SolveWorkspace {
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor|Eigen::AutoAlign> Matrix;
    Eigen::SVD<Matrix> svd;
    Eigen::VectorXd scratch;
};

struct BallCoefficientPrior::IntegrateWorkspace {
    int n1;
    int n2;
    double radiusSquared;
    double volumeFactor;
    Eigen::MatrixXd q;
    Eigen::VectorXd s;
    Eigen::VectorXd mu;
    Eigen::VectorXd nu;

    Eigen::VectorXd qtc;

    typedef Eigen::BlockReturnType<Eigen::MatrixXd>::Type MatrixBlock;
    typedef Eigen::BlockReturnType<Eigen::VectorXd>::SubVectorType VectorBlock;

    MatrixBlock getQ1() { return q.block(0, 0, q.rows(), n1); }
    MatrixBlock getQ2() { return q.block(0, n1, q.rows(), n2); }
    VectorBlock getQTC1() { return qtc.segment(0, n1); }
    VectorBlock getQTC2() { return qtc.segment(n1, n2); }
    VectorBlock getS1() { return s.segment(0, n1); }
    VectorBlock getS2() { return s.segment(n1, n2); }
    VectorBlock getNu1() { return nu.segment(0, n1); }
    VectorBlock getNu2() { return nu.segment(n1, n2); }

    bool sampleGaussian(Random & engine, double & weight) {
        drawDiagonalGaussian(engine, qtc, s, nu);
        if (qtc.squaredNorm() > radiusSquared) {
            weight = 0.0;
            return false;
        }
        return true;
    }

    bool sampleHybrid(Random & engine, double & weight) {
        drawDiagonalGaussian(engine, getQTC1(), getS1(), getNu1());
        double r1Squared = getQTC1().squaredNorm();
        if (r1Squared > radiusSquared) {
            weight = 0.0;
            return false;
        }
        double r2 = std::sqrt(radiusSquared - r1Squared);
        drawUniformOnSphere(engine, getQTC2(), r2);
        weight *= volumeFactor * std::pow(r2, 0.5 * n2);
        return true;
    }

};

double BallCoefficientPrior::operator()(
    ndarray::Array<Pixel const,1,1> const & coefficients
) const {
    if (ndarray::viewAsEigen(coefficients).norm() < _radius)
        return _normalization;
    return 0.0;
}

double BallCoefficientPrior::integrate(
    Random & engine,
    ndarray::Array<Pixel,2,2> const & coefficients,
    ndarray::Array<Pixel,1,1> const & weights,
    ndarray::Array<double const,2,2> const & modelMatrix,
    ndarray::Array<double const,1,1> const & dataVector
) const {
    typedef bool (IntegrateWorkspace::*SampleMethod)(Random &, double &);
    static int const MAX_REPLACEMENT_FACTOR = 2;
    ndarray::EigenView<Pixel,2,2> coeff(coefficients);
    if (!_integrateWorkspace) {
        _integrateWorkspace.reset(new IntegrateWorkspace());
    }
    IntegrateWorkspace & ws = *_integrateWorkspace; // just a short alias for legibility
    double constant = std::exp(-_solve(ws.mu, ws.q, ws.s, ws.n1, modelMatrix, dataVector));
    ws.radiusSquared = _radius * _radius;
    ws.qtc.resize(ws.mu.size());
    ws.n2 = ws.q.cols() - ws.n1;
    ws.nu = ws.q.transpose() * ws.mu;
    weights.deep() = 1.0;
    SampleMethod method;
    if (ws.n2 == 0) {
        // Likelihood Gaussian is nondegenerate; we sample from the likelihood.
        method = &IntegrateWorkspace::sampleGaussian;
    } else {
        // Likelihood Gaussian is degenerate; we draw from a hybrid Gaussian/uniform distribution.
        ws.volumeFactor = computeVolumeFactor(ws.n2);
        method = &IntegrateWorkspace::sampleHybrid;
    }
    int const mSamples = weights.getSize<0>();
    int mValid = 0;
    int mTotal = 0;
    while (mTotal < mSamples) {
        if ((ws.*method)(engine, weights[mValid])) {
            ndarray::viewAsEigen(coefficients[mValid]) = ws.q * ws.qtc;
            weights[mValid] *= _evaluate(coefficients[mValid]);
            if (weights[mValid] > std::numeric_limits<double>::epsilon()) ++mValid;
        }
        ++mTotal;
    }

    // Because (I think) sampling here should be relatively inexpensive and we're constrained
    // by how many samples we can put in a table, we should keep sampling until all samples
    // have nonzero weight, and just track how many we've rejected with mTotal.
    if (mValid * MAX_REPLACEMENT_FACTOR >= mSamples || mValid + 2 >= mSamples) {
        // Enough of our points are landing on nonzero prior points; we'll just loop through
        // and replace them with new Gaussian/Hybrid samples until we have a full set.
        while (mValid < mSamples) {
            if ((ws.*method)(engine, weights[mValid])) {
                ndarray::viewAsEigen(coefficients[mValid]) = ws.q * ws.qtc;
                weights[mValid] *= _evaluate(coefficients[mValid]);
                if (weights[mValid] > std::numeric_limits<double>::epsilon()) ++mValid;
            }
            ++mTotal;
        }
        assert(mValid == mSamples);
        double result = ndarray::viewAsEigen(weights).sum();
        ndarray::viewAsEigen(weights) /= result;
        result *= std::exp(0.5 * ws.n1 * std::log(2.0 * M_PI) - ws.getS1().cwise().log().sum())
            * constant / mTotal;
        return result;
    }

    // Most of our points are landing on points where the prior is zero; we'll sample from a
    // uniform distribution on a sphere and evaluate the Gaussian likelihood to fill in the rest.
    // We'll then do a variance-weighted mean of the Gaussian/Hybrid sample and the uniform sample.
    // In this case, we'll accept a point with zero weight, because sometimes the actual prior
    // may get much smaller than the spherical bound, and that means a noisy estimate is
    // better than an extremely slow one (it's a region of low probability anyway).
    if (mValid <= 2) {  // So small we can't even do variance weighting - just use uniform samples only
        mValid = 0;
    }
    int mUniform = mSamples - mValid;
    double volume = computeVolume(coefficients.getSize<1>(), _radius);
    for (int m = mValid; m < mSamples; ++m) {
        drawUniformOnSphere(engine, ws.qtc, _radius);
        ndarray::viewAsEigen(coefficients[m]) = ws.q * ws.qtc;
        ws.qtc -= ws.nu;
        ws.qtc.cwise() *= ws.s;
        weights[m] *= _evaluate(coefficients[m]) * volume * std::exp(-0.5 * ws.qtc.squaredNorm());
    }
    ndarray::EigenView<Pixel,1,1> w(weights);
    if (mUniform == mSamples) {
        double result = w.sum();
        w /= result;
        result /= mUniform;
        return result;
    } else {
        double gMean = w.segment(0, mValid).sum() / mTotal;
        double gVar = (w.segment(0, mValid).cwise() - gMean).cwise().square().sum() 
            / (mTotal * (mTotal - 1.0));
        double uMean = w.segment(mValid, mUniform).sum() / mTotal;
        double uVar = (w.segment(mValid, mUniform).cwise() - uMean).cwise().square().sum() 
            / (mUniform * (mUniform - 1.0));
        double var = 1.0 / (1.0 / gVar + 1.0 / uVar);
        double gf = var / gVar;
        double uf = var / uVar;
        w.segment(0, mValid) *= gf;
        w.segment(mValid, mUniform) *= uf;
        w /= w.sum();
        return constant * (gf * gMean + uf * uMean);
    }
}

BallCoefficientPrior::~BallCoefficientPrior() {}

double BallCoefficientPrior::computeVolume(int d, double radius, double factor) {
    if (factor <= 0.0) {
        factor = computeVolumeFactor(d);
    }
    return factor * std::pow(radius, d);
}

double BallCoefficientPrior::computeVolumeFactor(int d) {
    return std::pow(M_PI, 0.5 * d) / boost::math::tgamma(0.5 * d + 1.0);
}

double BallCoefficientPrior::_solve(
    Eigen::VectorXd & mu, Eigen::MatrixXd & q, Eigen::VectorXd & s, int & n1,
    ndarray::Array<double const,2,2> const & modelMatrix,
    ndarray::Array<double const,1,1> const & dataVector        
) const {
    static double const rcond = std::sqrt(std::numeric_limits<double>::epsilon());
    if (!_solveWorkspace) {
        _solveWorkspace.reset(new SolveWorkspace());
    }
    SolveWorkspace & ws = *_solveWorkspace; // just a short alias for legibility
    ws.svd.compute(ndarray::viewAsEigen(modelMatrix));
    ws.svd.sort();
    s = ws.svd.singularValues();
    q = ws.svd.matrixV();
    ws.scratch.resize(s.size());
    double sMin = s[0] * rcond;
    for (n1 = 1; n1 < s.size(); ++n1) {
        if (s.coeff(n1) <= sMin) {
            break;
        }
    }
    ws.scratch.segment(0, n1) = ws.svd.matrixU().block(0, 0, ws.svd.matrixU().rows(), n1).transpose() 
        * ndarray::viewAsEigen(dataVector);
    double result = 0.5 * (
        ws.svd.matrixU().block(0, 0, ws.svd.matrixU().rows(), n1) * ws.scratch.segment(0, n1)
        - ndarray::viewAsEigen(dataVector)
    ).squaredNorm();
    ws.scratch.segment(0, n1).cwise() /= s.segment(0, n1);
    mu = q.block(0, 0, q.rows(), n1) * ws.scratch;
    return result;
}

}}} // namespace lsst::meas::multifit
