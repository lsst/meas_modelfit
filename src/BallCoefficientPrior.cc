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


class BallCoefficientPrior::Impl {
public:
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor|Eigen::AutoAlign> Matrix;
    typedef Eigen::VectorXd Vector;
    typedef Eigen::BlockReturnType<Eigen::VectorXd>::SubVectorType SubVector;
    typedef Eigen::MatrixXd UMatrix;
    typedef Eigen::MatrixXd VMatrix;
    typedef Eigen::BlockReturnType<UMatrix>::Type UBlock;
    typedef Eigen::BlockReturnType<VMatrix>::Type VBlock;
    typedef Eigen::BlockReturnType< ndarray::EigenView<Pixel,1,1> >::SubVectorType CoeffView;

    typedef bool (Impl::*SampleMethod)(
        Random & engine, double & weight, ndarray::Array<Pixel,1,1> const & coeff
    );

    explicit Impl(BallCoefficientPrior const & parent) : _constraints(parent._constraints) {}

    int getN() const { return _s.size(); }

    int getN_a() const { return _n_a; }
    int getN_b() const { return _n_b; }

    Vector & getS() { return _s; }
    SubVector getS_a() { return _s.segment(0, _n_a); }
    SubVector getS_b() { return _s.segment(_n_a, _n_b); }

    Vector & getNu() { return _nu; }
    SubVector getNu_a() { return _nu.segment(0, _n_a); }
    SubVector getNu_b() { return _nu.segment(_n_a, _n_b); }

    Vector & getMu() { return _mu; }
    SubVector getMu(Constraint const & c) { return _mu.segment(c.offset, c.size); }

    UMatrix const & getU() { return _svd.matrixU(); }
    UBlock const getU_a() { return getU().block(0, 0, getU().rows(), _n_a); }
    UBlock const getU_b() { return getU().block(0, _n_a, getU().rows(), _n_b); }

    VMatrix const & getV() { return _svd.matrixV(); }
    VBlock const getV_a() { return getV().block(0, 0, getV().rows(), _n_a); }
    VBlock const getV_b() { return getV().block(0, _n_a, getV().rows(), _n_b); }
    VBlock const getV_a(Constraint const & c) { return getV().block(c.offset, 0, c.size, _n_a); }
    VBlock const getV_b(Constraint const & c) { return getV().block(c.offset, _n_a, c.size, _n_b); }

    Vector & getScratch() { return _scratch; }
    SubVector getScratch_a() { return _scratch.segment(0, _n_a); }
    SubVector getScratch_b() { return _scratch.segment(_n_a, _n_b); }

    bool sampleGaussian(Random & engine, double & weight, ndarray::Array<Pixel,1,1> const & coeff) {
        drawDiagonalGaussian(engine, getScratch(), getS(), getNu());
        ndarray::viewAsEigen(coeff) = getV() * getScratch();
        for (ConstraintIter i = _constraints.begin(); i != _constraints.end(); ++i) {
            if (ndarray::viewAsEigen(coeff).segment(i->offset, i->size).norm() > i->radius) {
                weight = 0.0;
                return false;
            }
        }
        return true;
    }

    bool sampleHybrid(Random & engine, double & weight, ndarray::Array<Pixel,1,1> const & coeff) {
        drawDiagonalGaussian(engine, getScratch_a(), getS_a(), getNu_a());
        ndarray::viewAsEigen(coeff) = getV_a() * getScratch_a();
        double r2Sum = 0.0;
        for (ConstraintIter i = _constraints.begin(); i != _constraints.end(); ++i) {
            double r2 = i->radius * i->radius 
                - ndarray::viewAsEigen(coeff).segment(i->offset, i->size).squaredNorm();
            if (r2 > 0.0) {
                weight = 0.0;
                return false;
            }
            r2Sum += r2;
        }
        double rSum = std::sqrt(r2Sum);
        drawUniformOnSphere(engine, getScratch_b(), rSum);
        ndarray::viewAsEigen(coeff) += getV_b() * getScratch_b();
        for (ConstraintIter i = _constraints.begin(); i != _constraints.end(); ++i) {
            if (ndarray::viewAsEigen(coeff).segment(i->offset, i->size).norm() > i->radius) {
                weight = 0.0;
                return false;
            }
        }
        weight *= computeVolume(getN_b(), rSum);
        return true;
    }

    void sampleConstraint(Random & engine, double & weight, ndarray::Array<Pixel,1,1> const & coeff) {
        for (ConstraintIter i = _constraints.begin(); i != _constraints.end(); ++i) {
            drawUniformOnSphere(engine, ndarray::viewAsEigen(coeff).segment(i->offset, i->size), i->radius);
        }
        getScratch_a() = getV_a().transpose() * ndarray::viewAsEigen(coeff);
        getScratch_a() -= getNu_a();
        getScratch_a() *= getS_a();
        weight *= std::exp(-0.5 * getScratch_a().squaredNorm());
    }

    double solve(
        ndarray::Array<double const,2,2> const & modelMatrix,
        ndarray::Array<double const,1,1> const & dataVector        
    );

private:
    ConstraintList const & _constraints;
    int _n_a;
    int _n_b;
    Vector _s;
    Vector _mu;
    Vector _nu;
    Vector _scratch;
    Eigen::SVD<Matrix> _svd;
};

double BallCoefficientPrior::Impl::solve(
    ndarray::Array<double const,2,2> const & modelMatrix,
    ndarray::Array<double const,1,1> const & dataVector        
) {
    static double const EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    _svd.compute(ndarray::viewAsEigen(modelMatrix));
    _svd.sort();
    _s = _svd.singularValues();
    _scratch.resize(_s.size());
    double sMin = _s[0] * EPS;
    for (_n_a = 1; _n_a < _s.size(); ++_n_a) {
        if (_s.coeff(_n_a) <= sMin) {
            break;
        }
    }
    _n_b = getN() - _n_a;
    getNu_a() = getU_a().transpose() * ndarray::viewAsEigen(dataVector);
    double result = 0.5 * (getU_a() * getNu_a() - ndarray::viewAsEigen(dataVector)).squaredNorm();
    getNu_a().cwise() /= getS_a();
    getNu_b().setZero();
    getMu() = getV_a() * getNu_a();
    return result;
}

double BallCoefficientPrior::operator()(lsst::ndarray::Array<Pixel const,1,1> const & coefficients) const {
    ndarray::EigenView<Pixel const,1,1> coeff(coefficients);
    for (ConstraintIter i = _constraints.begin(); i != _constraints.end(); ++i) {
        if (coeff.segment(i->offset, i->size).norm() > i->radius) {
            return 0.0;
        }
    }
    return evaluate(coefficients);
}

double BallCoefficientPrior::integrate(
    Random & engine,
    ndarray::Array<Pixel,2,2> const & coefficients,
    ndarray::Array<Pixel,1,1> const & weights,
    ndarray::Array<double const,2,2> const & modelMatrix,
    ndarray::Array<double const,1,1> const & dataVector
) const {
    static int const MAX_REPLACEMENT_FACTOR = 2;
    ndarray::EigenView<Pixel,2,2> coeff(coefficients);
    double constant = std::exp(-_impl->solve(modelMatrix, dataVector));
    weights.deep() = 1.0;
    Impl::SampleMethod method;
    if (_impl->getN_b() == 0) {
        // Likelihood Gaussian is nondegenerate; we sample from the likelihood.
        method = &Impl::sampleGaussian;
    } else {
        // Likelihood Gaussian is degenerate; we draw from a hybrid Gaussian/uniform distribution.
        method = &Impl::sampleHybrid;
    }
    double const likelihoodSampleVolume 
        = std::exp(0.5 * _impl->getN_a() * std::log(2.0 * M_PI) - _impl->getS_a().cwise().log().sum());
    int const mSamples = weights.getSize<0>();
    int mValid = 0;
    int mTotal = 0;
    while (mTotal < mSamples) {
        if (((*_impl).*method)(engine, weights[mValid], coefficients[mValid])) {
            weights[mValid] *= evaluate(coefficients[mValid]);
            if (weights[mValid] > std::numeric_limits<double>::epsilon()) ++mValid;
        }
        ++mTotal;
    }

    // Because sampling here should be relatively inexpensive and
    // we're constrained by how many samples we can put in a table, we
    // should keep sampling until all samples have nonzero weight, and
    // just track how many we've rejected with mTotal.
    if (mValid * MAX_REPLACEMENT_FACTOR >= mSamples || mValid + 2 >= mSamples) {
        // Enough of our points are landing on nonzero prior points;
        // we'll just loop through and replace them with new
        // Gaussian/Hybrid samples until we have a full set.
        while (mValid < mSamples) {
            if (((*_impl).*method)(engine, weights[mValid], coefficients[mValid])) {
                weights[mValid] *= evaluate(coefficients[mValid]);
                if (weights[mValid] > std::numeric_limits<double>::epsilon()) ++mValid;
            }
            ++mTotal;
        }
        assert(mValid == mSamples);
        double result = ndarray::viewAsEigen(weights).sum();
        ndarray::viewAsEigen(weights) /= result;
        result *= constant * likelihoodSampleVolume / mTotal;
        return result;
    }

    // Most of our points are landing on points where the prior is
    // zero; we'll fill the rest of the table by sampling from the prior and evaluating the likelihood,
    // then do an inverse-variance weighted mean of the two methods.
    ndarray::EigenView<Pixel,1,1> w(weights);
    if (mValid <= 2) {  // Too small to do variance weighting - just use prior samples only
        mValid = 0;
    } else {
        w.segment(0, mValid) *= likelihoodSampleVolume;
    }
    int mUniform = mSamples - mValid;
    for (int m = mValid; m < mSamples; ++m) {
        _impl->sampleConstraint(engine, weights[m], coefficients[m]);
        weights[m] *= evaluate(coefficients[m]);
    }
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

double BallCoefficientPrior::computeVolume(int d, double radius) {
    return computeVolumeFactor(d) * std::pow(radius, d);
}

double BallCoefficientPrior::computeVolumeFactor(int d) {
    static int const N = 20;
    static double const cached[N] = {
        1.0,
        2.0,
        3.1415926535897931,
        4.1887902047863914,
        4.934802200544679,
        5.263789013914324,
        5.1677127800499694,
        4.7247659703314007,
        4.0587121264167676,
        3.2985089027387064,
        2.5501640398773451,
        1.8841038793898999,
        1.3352627688545893,
        0.91062875478328287,
        0.59926452932079188,
        0.38144328082330436,
        0.23533063035889312,
        0.140981106917139,
        0.082145886611128191,
        0.046621601030088528
    };
    return (d < N) ? cached[d] : std::pow(M_PI, 0.5 * d) / boost::math::tgamma(0.5 * d + 1.0);
}

}}} // namespace lsst::meas::multifit
