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

#include "boost/math/special_functions/gamma.hpp"

#include "ndarray/eigen.h"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/Persistable.cc"
#include "lsst/meas/modelfit/Mixture.h"

namespace tbl = lsst::afw::table;

namespace lsst {
namespace afw {
namespace table {
namespace io {

template std::shared_ptr<meas::modelfit::Mixture> PersistableFacade<meas::modelfit::Mixture>::dynamicCast(
        std::shared_ptr<Persistable> const&);

}  // namespace io
}  // namespace table
}  // namespace afw
namespace meas {
namespace modelfit {

void MixtureComponent::setSigma(Matrix const & sigma) {
    _sigmaLLT.compute(sigma);
    _sqrtDet = _sigmaLLT.matrixLLT().diagonal().prod();
}

MixtureComponent MixtureComponent::project(int dim) const {
    Vector mu(1);
    mu << _mu[dim];
    Matrix sigma(1,1);
    sigma << getSigma()(dim, dim);
    return MixtureComponent(weight, mu, sigma);
}

MixtureComponent MixtureComponent::project(int dim1, int dim2) const {
    Vector mu(2);
    mu << _mu[dim1], _mu[dim2];
    Matrix sigma(2,2);
    Matrix fullSigma = getSigma();
    sigma <<
        fullSigma(dim1, dim1), fullSigma(dim1, dim2),
        fullSigma(dim2, dim1), fullSigma(dim2, dim2);
    return MixtureComponent(weight, mu, sigma);
}

MixtureComponent::MixtureComponent(int dim) :
    weight(1.0), _sqrtDet(1.0), _mu(Vector::Zero(dim)), _sigmaLLT(Matrix::Identity(dim,dim)) {}


MixtureComponent::MixtureComponent(Scalar weight_, Vector const & mu, Matrix const & sigma) :
    weight(weight_), _mu(mu), _sigmaLLT(mu.size())
{
    LSST_THROW_IF_NE(
        sigma.rows(), _mu.size(),
        pex::exceptions::LengthError,
        "Number of rows of sigma matrix (%d) does not match size of mu vector (%d)"
    );
    LSST_THROW_IF_NE(
        sigma.cols(), _mu.size(),
        pex::exceptions::LengthError,
        "Number of columns of sigma matrix (%d) does not match size of mu vector (%d)"
    );
    _sigmaLLT.compute(sigma);
    _sqrtDet = _sigmaLLT.matrixLLT().diagonal().prod();
}

MixtureComponent & MixtureComponent::operator=(MixtureComponent const & other) {
    LSST_THROW_IF_NE(
        other.getDimension(), getDimension(),
        pex::exceptions::LengthError,
        "Cannot assign MixtureComponent with dim=%d to one with dim=%d"
    );
    if (&other != this) {
        _sqrtDet = other._sqrtDet;
        _mu = other._mu;
        _sigmaLLT = other._sigmaLLT;
    }
    return *this;
}

void MixtureComponent::_stream(std::ostream & os, int offset) const {
    static Eigen::IOFormat muFormat(12, 0, ",", "\n", "[", "]", "[", "]");
    std::string pad(offset, ' ');
    Eigen::IOFormat sigmaFormat(12, 0, ",", ",\n        " + pad, "[", "]", "[", "]");
    os << pad << "( weight=" << weight << ", mu="
       << _mu.transpose().format(muFormat) << "," << std::endl;
    os << pad << "  sigma=" << getSigma().format(sigmaFormat) << " )" << std::endl;
}

std::shared_ptr<Mixture> Mixture::project(int dim) const {
    ComponentList components;
    components.reserve(size());
    for (const_iterator i = begin(); i != end(); ++i) {
        components.push_back(i->project(dim));
    }
    return std::make_shared<Mixture>(1, components, _df);
}

std::shared_ptr<Mixture> Mixture::project(int dim1, int dim2) const {
    ComponentList components;
    components.reserve(size());
    for (const_iterator i = begin(); i != end(); ++i) {
        components.push_back(i->project(dim1, dim2));
    }
    return std::make_shared<Mixture>(2, components, _df);
}

void Mixture::normalize() {
    Scalar sum = 0.0;
    for (iterator i = begin(); i != end(); ++i) {
        sum += i->weight;
    }
    for (iterator i = begin(); i != end(); ++i) {
        i->weight /= sum;
    }
}

void Mixture::shift(int dim, Scalar offset) {
    for (iterator i = begin(); i != end(); ++i) {
        i->_mu[dim] += offset;
    }
}

std::size_t Mixture::clip(Scalar threshold) {
    std::size_t count = 0;
    iterator i = begin();
    while (i != end()) {
        if (i->weight <= threshold) {
            i = _components.erase(i);
            ++count;
        } else {
            ++i;
        }
    }
    if (count) normalize();
    return count;
}

void Mixture::setDegreesOfFreedom(Scalar df) {
    _df = df;
    if (_df == std::numeric_limits<Scalar>::infinity()) {
        _norm = std::pow(2.0 * M_PI, 0.5*_dim);
        _isGaussian = true;
    } else {
        _norm = boost::math::tgamma_delta_ratio(0.5*_df, 0.5*_dim) * std::pow(_df*M_PI, 0.5*_dim);
        _isGaussian = false;
    }
}

void Mixture::evaluate(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar,1,0> const & p
) const {
    LSST_THROW_IF_NE(
        x.getSize<0>(), p.getSize<0>(),
        pex::exceptions::LengthError,
        "First dimension of x array (%d) does not match size of p array (%d)"
    );
    LSST_THROW_IF_NE(
        x.getSize<1>(), _dim,
        pex::exceptions::LengthError,
        "Second dimension of x array (%d) does not dimension of mixture (%d)"
    );
    ndarray::Array<Scalar const,2,1>::Iterator ix = x.begin(), xEnd = x.end();
    ndarray::Array<Scalar,1,0>::Iterator ip = p.begin();
    for (; ix != xEnd; ++ix, ++ip) {
        *ip = evaluate(ndarray::asEigenMatrix(*ix));
    }
}

void Mixture::evaluateComponents(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar,2,1> const & p
) const {
    LSST_THROW_IF_NE(
        x.getSize<0>(), p.getSize<0>(),
        pex::exceptions::LengthError,
        "First dimension of x array (%d) does not match first dimension of p array (%d)"
    );
    LSST_THROW_IF_NE(
        x.getSize<1>(), _dim,
        pex::exceptions::LengthError,
        "Second dimension of x array (%d) does not dimension of mixture (%d)"
    );
    LSST_THROW_IF_NE(
        p.getSize<1>(), static_cast<int>(_components.size()),
        pex::exceptions::LengthError,
        "Second dimension of p array (%d) does not match number of components (%d)"
    );
    ndarray::Array<Scalar const,2,1>::Iterator ix = x.begin(), xEnd = x.end();
    ndarray::Array<Scalar,2,1>::Iterator ip = p.begin();
    for (; ix != xEnd; ++ix, ++ip) {
        ndarray::Array<Scalar,2,1>::Reference::Iterator jp = ip->begin();
        for (const_iterator j = begin(); j != end(); ++j, ++jp) {
            *jp = evaluate(*j, ndarray::asEigenMatrix(*ix));
        }
    }
}

void Mixture::evaluateDerivatives(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> & x,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> & gradient,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> & hessian
) const {
    _evaluateDerivativesImpl(x, gradient, &hessian);
}

void Mixture::evaluateDerivatives(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> & x,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> & gradient
) const {
    _evaluateDerivativesImpl<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>,
                             Eigen::Matrix<Scalar, Eigen::Dynamic, 1>,
                             Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> (x, gradient, nullptr, false);
}

void Mixture::evaluateDerivatives(
    ndarray::Array<Scalar const,1,1> const & x,
    ndarray::Array<Scalar,1,1> const & gradient,
    ndarray::Array<Scalar,2,1> const & hessian
) const {
    auto hessianEigen = ndarray::asEigenMatrix(hessian);
    auto gradientEigen = ndarray::asEigenMatrix(gradient);
    _evaluateDerivativesImpl(ndarray::asEigenMatrix(x), gradientEigen, &hessianEigen);
}

template <typename A, typename B, typename C>
void Mixture::_evaluateDerivativesImpl(
    A const & x,
    B & gradient,
    C * hessian,
    bool computeHessian
) const {
    LSST_THROW_IF_NE(
        x.size(), _dim,
        pex::exceptions::LengthError,
        "Size of x array (%d) does not dimension of mixture (%d)"
    );
    LSST_THROW_IF_NE(
        gradient.rows(), _dim,
        pex::exceptions::LengthError,
        "Size of gradient array (%d) does not dimension of mixture (%d)"
    );
    gradient.setZero();

    if (computeHessian) {
        if (!hessian) {
            LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                "Pointer to hessian object must be specified if computeHessian is true"
           );
        }
        size_t rows = hessian->rows();
        size_t columns = hessian->cols();
        LSST_THROW_IF_NE(
            rows, _dim,
            pex::exceptions::LengthError,
            "Number of rows of hessian array (%d) does not dimension of mixture (%d)"
        );
        LSST_THROW_IF_NE(
            columns, _dim,
            pex::exceptions::LengthError,
            "Number of columns of hessian array (%d) does not dimension of mixture (%d)"
        );
        hessian->setZero();
    }
    Eigen::MatrixXd sigmaInv(_dim, _dim);
    for (ComponentList::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        _workspace = x - i->_mu;
        i->_sigmaLLT.matrixL().solveInPlace(_workspace);
        Scalar z = _workspace.squaredNorm();
        i->_sigmaLLT.matrixL().adjoint().solveInPlace(_workspace);
        sigmaInv.setIdentity();
        i->_sigmaLLT.matrixL().solveInPlace(sigmaInv);
        i->_sigmaLLT.matrixL().adjoint().solveInPlace(sigmaInv);
        Scalar f = _evaluate(z) / i->_sqrtDet;
        if (_isGaussian) {
            gradient += -i->weight * f * _workspace;
            if (computeHessian) {
                *hessian += i->weight * f * (_workspace * _workspace.adjoint() - sigmaInv);
            }
        } else {
            double v = (_dim + _df) / (_df + z);
            double u = v*v*(1.0 + 2.0/(_dim + _df));
            gradient += -i->weight * f * v * _workspace;
            if (computeHessian) {
                *hessian += i->weight * f * (u * _workspace * _workspace.adjoint() - v * sigmaInv);
            }
        }
    }
}

void Mixture::draw(afw::math::Random & rng, ndarray::Array<Scalar,2,1> const & x) const {
    ndarray::Array<Scalar,2,1>::Iterator ix = x.begin(), xEnd = x.end();
    std::vector<Scalar> cumulative;
    cumulative.reserve(_components.size());
    Scalar sum = 0.0;
    for (const_iterator k = begin(); k != end(); ++k) {
        sum += k->weight;
        cumulative.push_back(sum);
    }
    cumulative.back() = 1.0;
    for (; ix != xEnd; ++ix) {
        Scalar target = rng.uniform();
        std::size_t k = std::lower_bound(cumulative.begin(), cumulative.end(), target)
            - cumulative.begin();
        assert(k != cumulative.size());
        Component const & component = _components[k];
        for (int j = 0; j < _dim; ++j) {
            _workspace[j] = rng.gaussian();
        }
        if (!_isGaussian) {
            _workspace *= std::sqrt(_df/rng.chisq(_df));
        }
        ndarray::asEigenMatrix(*ix) = component._mu + (component._sigmaLLT.matrixL() * _workspace);
    }
}

void Mixture::updateEM(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar const,1,0> const & w,
    UpdateRestriction const & restriction,
    Scalar tau1, Scalar tau2
) {
    LSST_THROW_IF_NE(
        x.getSize<0>(), w.getSize<0>(),
        pex::exceptions::LengthError,
        "First dimension of x array (%d) does not match size of w array (%d)"
    );
    LSST_THROW_IF_NE(
        x.getSize<1>(), _dim,
        pex::exceptions::LengthError,
        "Second dimension of x array (%d) does not dimension of mixture (%d)"
    );
    int const nSamples = w.getSize<0>();
    int const nComponents = _components.size();
    Matrix p(nSamples, nComponents);
    Matrix gamma(nSamples, nComponents);
    for (int i = 0; i < nSamples; ++i) {
        Scalar pSum = 0.0;
        for (int k = 0; k < nComponents; ++k) {
            double z = _computeZ(_components[k], ndarray::asEigenMatrix(x[i]));
            pSum += p(i, k) = _components[k].weight*_evaluate(z)/_components[k]._sqrtDet;
            if (!_isGaussian) {
                gamma(i, k) = (_df + _dim) / (_df + z);
            }
        }
        p.row(i) *= w[i] / pSum;
    }
    if (_isGaussian) {
        for (int k = 0; k < nComponents; ++k) {
            double weight = _components[k].weight = p.col(k).sum();
            Vector & mu = _components[k]._mu;
            Matrix sigma = Matrix::Zero(_dim, _dim);
            mu = (p.col(k).adjoint() * ndarray::asEigenMatrix(x)) / weight;
            restriction.restrictMu(mu);
            Vector dx = Vector::Zero(_dim);
            for (int i = 0; i < nSamples; ++i) {
                dx = ndarray::asEigenMatrix(x[i]) - mu;
                sigma.selfadjointView<Eigen::Lower>().rankUpdate(dx, p(i, k));
            }
            sigma /= weight;
            restriction.restrictSigma(sigma);
            updateDampedSigma(k, sigma, tau1, tau2);
        }
    } else {
        for (int k = 0; k < nComponents; ++k) {
            double weight = _components[k].weight = p.col(k).sum();
            Vector & mu = _components[k]._mu;
            Matrix sigma = Matrix::Zero(_dim, _dim);
            mu =
                ((p.col(k).array() * gamma.col(k).array()).matrix().adjoint() * ndarray::asEigenMatrix(x))
                / p.col(k).dot(gamma.col(k));
            restriction.restrictMu(mu);
            Vector dx = Vector::Zero(_dim);
            for (int i = 0; i < nSamples; ++i) {
                dx = ndarray::asEigenMatrix(x[i]) - mu;
                sigma.selfadjointView<Eigen::Lower>().rankUpdate(dx, gamma(i, k) * p(i, k));
            }
            sigma /= weight;
            restriction.restrictSigma(sigma);
            updateDampedSigma(k, sigma, tau1, tau2);
        }
    }
}

void Mixture::updateEM(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar const,1,0> const & w,
    Scalar tau1, Scalar tau2
) {
    updateEM(x, w, UpdateRestriction(_dim), tau1, tau2);
}

void Mixture::updateEM(
    ndarray::Array<Scalar const,2,1> const & x,
    UpdateRestriction const & restriction,
    Scalar tau1, Scalar tau2
) {
    ndarray::Array<Scalar,1,1> w = ndarray::allocate(x.getSize<0>());
    w.deep() = 1.0 / w.getSize<0>();
    updateEM(x, w, restriction, tau1, tau2);
}

std::shared_ptr<Mixture> Mixture::clone() const {
    return std::make_shared<Mixture>(*this);
}

Mixture::Mixture(int dim, ComponentList & components, Scalar df) :
    _dim(dim), _df(0.0), _workspace(dim)
{
    setDegreesOfFreedom(df);
    _components.swap(components);
    normalize();
}

void Mixture::updateDampedSigma(int k, Matrix const & sigma, double tau1, double tau2) {
    Eigen::LLT<Matrix> sigmaLLT(sigma);
    Scalar sqrtDet = sigmaLLT.matrixLLT().diagonal().prod();
    Scalar r = sqrtDet / _components[k]._sqrtDet;
    if (!(r >= tau1)) {
        Scalar beta1 = 2*(tau2 - 1.0)/tau1;
        Scalar beta2 = -0.5/tau1;
        Scalar alpha = beta1*r*(1.0 + beta2*r);
        _components[k].setSigma(alpha*sigma + (1.0 - alpha)*_components[k].getSigma());
    } else {
        _components[k]._sigmaLLT = sigmaLLT;
        _components[k]._sqrtDet = sqrtDet;
    }
}

Scalar Mixture::_evaluate(Scalar z) const {
    if (_isGaussian) {
        return std::exp(-0.5*z) / _norm;
    } else {
        return std::pow(z/_df + 1.0, -0.5*(_df + _dim)) / _norm;
    }
}

void Mixture::_stream(std::ostream & os) const {
    os << "Mixture(dim=" << _dim << ", [\n";
    for (const_iterator i = begin(); i != end(); ++i) {
        i->_stream(os, 2);
    }
    os << "],\n  df=" << _df << "\n)";
}

namespace {

class MixturePersistenceKeys {
public:
    tbl::Schema schema;
    tbl::Key<Scalar> weight;
    tbl::Key< tbl::Array<Scalar> > mu;
    tbl::Key< tbl::Array<Scalar> > sigma;

    explicit MixturePersistenceKeys(int dim) :
        schema(),
        weight(schema.addField<Scalar>("weight", "weight of mixture component")),
        mu(schema.addField<tbl::Array<Scalar> >("mu", "location parameter", dim)),
        sigma(schema.addField<tbl::Array<Scalar> >("sigma", "size/shape parameter", dim*dim))
    {}

    explicit MixturePersistenceKeys(tbl::Schema const & schema_) :
        schema(schema_),
        weight(schema["weight"]),
        mu(schema["mu"]),
        sigma(schema["sigma"])
    {}

    // No copying
    MixturePersistenceKeys (const MixturePersistenceKeys&) = delete;
    MixturePersistenceKeys& operator=(const MixturePersistenceKeys&) = delete;
    
    // No moving
    MixturePersistenceKeys (MixturePersistenceKeys&&) = delete;
    MixturePersistenceKeys& operator=(MixturePersistenceKeys&&) = delete;
};

class MixtureFactory : public tbl::io::PersistableFactory {
public:

    virtual std::shared_ptr<tbl::io::Persistable>
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() >= 1u);
        MixturePersistenceKeys const & keys = MixturePersistenceKeys(catalogs.front().getSchema());
        tbl::BaseCatalog::const_iterator iter = catalogs.front().begin();
        tbl::BaseCatalog::const_iterator const end = catalogs.front().end();
        int dim = keys.mu.getSize();
        double df = iter->get(keys.weight); // use first record to store df, nothing else.
        ++iter;
        Mixture::ComponentList components;
        components.reserve(end - iter);
        for (; iter != end; ++iter) {
            components.push_back(
                MixtureComponent(
                    iter->get(keys.weight),
                    Eigen::Map<Vector const>(iter->getElement(keys.mu), dim),
                    Eigen::Map<Matrix const>(iter->getElement(keys.sigma), dim, dim)
                )
            );
        }
        return std::make_shared<Mixture>(dim, components, df);
    }

    explicit MixtureFactory(std::string const & name) : tbl::io::PersistableFactory(name) {}

};

std::string getMixturePersistenceName() { return "Mixture"; }

MixtureFactory registration(getMixturePersistenceName());

} // anonymous

std::string Mixture::getPersistenceName() const { return getMixturePersistenceName(); }

void Mixture::write(OutputArchiveHandle & handle) const {
    MixturePersistenceKeys const keys(_dim);
    tbl::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    std::shared_ptr<tbl::BaseRecord> record = catalog.addNew();
    record->set(keys.weight, _df);
    for (const_iterator i = begin(); i != end(); ++i) {
        record = catalog.addNew();
        record->set(keys.weight, i->weight);
        Eigen::Map<Vector>(record->getElement(keys.mu), _dim) = i->_mu;
        Eigen::Map<Matrix>(record->getElement(keys.sigma), _dim, _dim) = i->getSigma();
    }
    handle.saveCatalog(catalog);
}

}}} // namespace lsst::meas::modelfit
