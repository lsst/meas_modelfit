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
#include "lsst/meas/multifit/Mixture.h"

namespace tbl = lsst::afw::table;

namespace lsst { namespace meas { namespace multifit {

template <int N>
void MixtureComponent<N>::_stream(std::ostream & os, int offset) const {
    static Eigen::IOFormat muFormat(12, 0, ",", "\n", "[", "]", "[", "]");
    std::string pad(offset, ' ');
    Eigen::IOFormat sigmaFormat(12, 0, ",", ",\n        " + pad, "[", "]", "[", "]");
    os << pad << "( weight=" << weight << ", mu="
       << _mu.transpose().format(muFormat) << "," << std::endl;
    os << pad << "  sigma=" << getSigma().format(sigmaFormat) << " )" << std::endl;
}

template <int N>
MixtureComponent<1> MixtureComponent<N>::project(int dim) const {
    MixtureComponent<1>::Vector mu;
    mu << _mu[dim];
    MixtureComponent<1>::Matrix sigma;
    sigma << getSigma()(dim, dim);
    return MixtureComponent<1>(weight, mu, sigma);
}

template <int N>
MixtureComponent<2> MixtureComponent<N>::project(int dim1, int dim2) const {
    MixtureComponent<2>::Vector mu;
    mu << _mu[dim1], _mu[dim2];
    MixtureComponent<2>::Matrix sigma;
    Matrix fullSigma = getSigma();
    sigma <<
        fullSigma(dim1, dim1), fullSigma(dim1, dim2),
        fullSigma(dim2, dim1), fullSigma(dim2, dim2);
    return MixtureComponent<2>(weight, mu, sigma);
}

template <int N>
void Mixture<N>::_stream(std::ostream & os) const {
    os << "Mixture([\n";
    for (const_iterator i = begin(); i != end(); ++i) {
        i->_stream(os, 2);
    }
    os << "],\n  df=" << _df << "\n)";
}

template <int N>
Mixture<1> Mixture<N>::project(int dim) const {
    Mixture<1>::ComponentList components;
    components.reserve(size());
    for (const_iterator i = begin(); i != end(); ++i) {
        components.push_back(i->project(dim));
    }
    return Mixture<1>(components, _df);
}

template <int N>
Mixture<2> Mixture<N>::project(int dim1, int dim2) const {
    Mixture<2>::ComponentList components;
    components.reserve(size());
    for (const_iterator i = begin(); i != end(); ++i) {
        components.push_back(i->project(dim1, dim2));
    }
    return Mixture<2>(components, _df);
}

template <int N>
void Mixture<N>::normalize() {
    Scalar sum = 0.0;
    for (iterator i = begin(); i != end(); ++i) {
        sum += i->weight;
    }
    for (iterator i = begin(); i != end(); ++i) {
        i->weight /= sum;
    }
}

template <int N>
std::size_t Mixture<N>::clip(Scalar threshold) {
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

template <int N>
void Mixture<N>::setDegreesOfFreedom(Scalar df) {
    _df = df;
    if (_df == std::numeric_limits<Scalar>::infinity()) {
        _norm = std::pow(2.0 * M_PI, 0.5*N);
        _isGaussian = true;
    } else {
        _norm = boost::math::tgamma_delta_ratio(0.5*_df, 0.5*N) * std::pow(_df*M_PI, 0.5*N);
        _isGaussian = false;
    }
}

template <int N>
void Mixture<N>::evaluate(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar,1,0> const & p
) const {
    LSST_ASSERT_EQUAL(
        x.getSize<0>(), p.getSize<0>(),
        "First dimension of x array (%d) does not match size of p array (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        x.getSize<1>(), N,
        "Second dimension of x array (%d) does not dimension of mixture (%d)",
        pex::exceptions::LengthErrorException
    );
    ndarray::Array<Scalar const,2,1>::Iterator ix = x.begin(), xEnd = x.end();
    ndarray::Array<Scalar,1,0>::Iterator ip = p.begin();
    for (; ix != xEnd; ++ix, ++ip) {
        *ip = evaluate(ix->asEigen<N,1>());
    }
}

template <int N>
void Mixture<N>::evaluateComponents(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar,2,1> const & p
) const {
    LSST_ASSERT_EQUAL(
        x.getSize<0>(), p.getSize<0>(),
        "First dimension of x array (%d) does not match first dimension of p array (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        x.getSize<1>(), N,
        "Second dimension of x array (%d) does not dimension of mixture (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        p.getSize<1>(), static_cast<int>(_components.size()),
        "Second dimension of p array (%d) does not match number of components (%d)",
        pex::exceptions::LengthErrorException
    );
    ndarray::Array<Scalar const,2,1>::Iterator ix = x.begin(), xEnd = x.end();
    ndarray::Array<Scalar,2,1>::Iterator ip = p.begin();
    for (; ix != xEnd; ++ix, ++ip) {
        ndarray::Array<Scalar,2,1>::Reference::Iterator jp = ip->begin();
        for (const_iterator j = begin(); j != end(); ++j, ++jp) {
            *jp = evaluate(*j, ix->asEigen<N,1>());
        }
    }
}

template <int N>
void Mixture<N>::draw(afw::math::Random & rng, ndarray::Array<Scalar,2,1> const & x) const {
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
        for (int j = 0; j < N; ++j) {
            _workspace[j] = rng.gaussian();
        }
        if (_isGaussian) {
            ix->asEigen<N,1>() = component._mu + (component._sigmaLLT.matrixL() * _workspace);
        } else {
            ix->asEigen<N,1>() = component._mu
                + std::sqrt(_df/rng.chisq(_df)) * (component._sigmaLLT.matrixL() * _workspace);
        }
    }
}

template <int N>
void Mixture<N>::updateEM(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar const,1,1> const & w,
    UpdateRestriction const & restriction,
    Scalar tau1, Scalar tau2
) {
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> MatrixXs;
    LSST_ASSERT_EQUAL(
        x.getSize<0>(), w.getSize<0>(),
        "First dimension of x array (%d) does not match size of w array (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        x.getSize<1>(), N,
        "Second dimension of x array (%d) does not dimension of mixture (%d)",
        pex::exceptions::LengthErrorException
    );
    int const nSamples = w.getSize<0>();
    int const nComponents = _components.size();
    MatrixXs p(nSamples, nComponents);
    MatrixXs gamma(nSamples, nComponents);
    for (int i = 0; i < nSamples; ++i) {
        Scalar pSum = 0.0;
        for (int k = 0; k < nComponents; ++k) {
            double z = _computeZ(_components[k], x[i].asEigen());
            pSum += p(i, k) = _components[k].weight*_evaluate(z)/_components[k]._sqrtDet;
            if (!_isGaussian) {
                gamma(i, k) = (_df + N) / (_df + z);
            }
        }
        p.row(i) *= w[i] / pSum;
    }
    if (_isGaussian) {
        for (int k = 0; k < nComponents; ++k) {
            double weight = _components[k].weight = p.col(k).sum();
            Vector & mu = _components[k]._mu;
            Matrix sigma = Matrix::Zero();
            mu = (p.col(k).adjoint() * x.asEigen()) / weight;
            restriction.restrictMu(mu);
            Vector dx = Vector::Zero();
            for (int i = 0; i < nSamples; ++i) {
                dx = x[i].asEigen() - mu;
                sigma.template selfadjointView<Eigen::Lower>().rankUpdate(dx, p(i, k));
            }
            sigma /= weight;
            restriction.restrictSigma(sigma);
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
    } else {
        for (int k = 0; k < nComponents; ++k) {
            double weight = _components[k].weight = p.col(k).sum();
            Vector & mu = _components[k]._mu;
            Matrix sigma = Matrix::Zero();
            mu =
                ((p.col(k).array() * gamma.col(k).array()).matrix().adjoint() * x.asEigen())
                / p.col(k).dot(gamma.col(k));
            restriction.restrictMu(mu);
            Vector dx = Vector::Zero();
            for (int i = 0; i < nSamples; ++i) {
                dx = x[i].asEigen() - mu;
                sigma.template selfadjointView<Eigen::Lower>().rankUpdate(dx, gamma(i, k) * p(i, k));
            }
            sigma /= weight;
            restriction.restrictSigma(sigma);
            _components[k].setSigma(sigma);
        }
    }
}

template <int N>
void Mixture<N>::updateEM(
    ndarray::Array<Scalar const,2,1> const & x,
    ndarray::Array<Scalar const,1,1> const & w,
    Scalar tau1, Scalar tau2
) {
    updateEM(x, w, UpdateRestriction(), tau1, tau2);
}

template <int N>
void Mixture<N>::updateEM(
    ndarray::Array<Scalar const,2,1> const & x,
    UpdateRestriction const & restriction,
    Scalar tau1, Scalar tau2
) {
    ndarray::Array<Scalar,1,1> w = ndarray::allocate(x.getSize<0>());
    w.deep() = 1.0 / w.getSize<0>();
    updateEM(x, w, restriction, tau1, tau2);
}

template <int N>
Mixture<N>::Mixture(ComponentList & components, Scalar df) :
    _df(0.0)
{
    setDegreesOfFreedom(df);
    _components.swap(components);
    normalize();
}

template <int N>
samples::Scalar Mixture<N>::_evaluate(Scalar z) const {
    if (_isGaussian) {
        return std::exp(-0.5*z) / _norm;
    } else {
        return std::pow(z/_df + 1.0, -0.5*(_df + N)) / _norm;
    }
}

template <int N>
struct MixturePersistenceName;

template <int N>
class MixturePersistenceKeys : private boost::noncopyable {
public:
    tbl::Schema schema;
    samples::ScalarKey weight;
    samples::ArrayKey mu;
    samples::ArrayKey sigma;

    static MixturePersistenceKeys const & get() {
        static MixturePersistenceKeys const instance;
        return instance;
    }

private:
    MixturePersistenceKeys() :
        schema(),
        weight(schema.addField<samples::Scalar>("weight", "weight of mixture component")),
        mu(schema.addField<samples::ArrayTag>("mu", "location parameter", N)),
        sigma(schema.addField<samples::ArrayTag>("sigma", "size/shape parameter", N*N))
    {
        schema.getCitizen().markPersistent();
    }
};

template <int N>
class MixtureFactory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        MixturePersistenceKeys<N> const & keys = MixturePersistenceKeys<N>::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() >= 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys.schema);
        tbl::BaseCatalog::const_iterator iter = catalogs.front().begin();
        tbl::BaseCatalog::const_iterator const end = catalogs.front().end();
        double df = iter->get(keys.weight); // use first record to store df, nothing else.
        ++iter;
        typename Mixture<N>::ComponentList components;
        components.reserve(end - iter);
        for (; iter != end; ++iter) {
            components.push_back(
                MixtureComponent<N>(
                    iter->get(keys.weight),
                    Eigen::Map<typename Mixture<N>::Vector const>(iter->getElement(keys.mu), N),
                    Eigen::Map<typename Mixture<N>::Matrix const>(iter->getElement(keys.sigma), N, N)
                )
            );
        }
        return PTR(Mixture<N>)(new Mixture<N>(components, df));
    }

    explicit MixtureFactory(std::string const & name) : tbl::io::PersistableFactory(name) {}

};

template <int N>
std::string Mixture<N>::getPersistenceName() const { return MixturePersistenceName<N>::get(); }

template <int N>
void Mixture<N>::write(OutputArchiveHandle & handle) const {
    MixturePersistenceKeys<N> const & keys = MixturePersistenceKeys<N>::get();
    tbl::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    PTR(tbl::BaseRecord) record = catalog.addNew();
    record->set(keys.weight, _df);
    for (const_iterator i = begin(); i != end(); ++i) {
        record = catalog.addNew();
        record->set(keys.weight, i->weight);
        Eigen::Map<Vector>(record->getElement(keys.mu), N) = i->_mu;
        Eigen::Map<Matrix>(record->getElement(keys.sigma), N, N) = i->getSigma();
    }
    handle.saveCatalog(catalog);
}

#define INSTANTIATE(N)                          \
    template <>                                 \
    struct MixturePersistenceName<N> {          \
        static std::string get() { return "Mixture" # N; }      \
    };                                                          \
    template class MixtureComponent<N>;         \
    template class MixtureUpdateRestriction<N>; \
    template class Mixture<N>;                  \
    template class MixturePersistenceKeys<N>;   \
    template class MixtureFactory<N>;           \
    static MixtureFactory<N> registration ## N(MixturePersistenceName<N>::get())

INSTANTIATE(1);
INSTANTIATE(2);
INSTANTIATE(3);

}}} // namespace lsst::meas::multifit
