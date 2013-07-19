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

#ifndef LSST_MEAS_MULTIFIT_Mixture_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Mixture_h_INCLUDED

#include <limits>

#include "Eigen/Cholesky"
#include "Eigen/StdVector"

#include "ndarray.h"

#include "lsst/base.h"
#include "lsst/afw/math/Random.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

template <int N> class MixtureComponent;
template <int N> class Mixture;

template <int N>
std::ostream & operator<<(std::ostream & os, MixtureComponent<N> const & self);

template <int N>
std::ostream & operator<<(std::ostream & os, Mixture<N> const & self);

/**
 *  @brief A weighted Student's T or Gaussian distribution used as a component in a Mixture.
 */
template <int N>
class MixtureComponent {
public:

    typedef samples::Scalar Scalar;
    typedef Eigen::Matrix<samples::Scalar,N,1> Vector;
    typedef Eigen::Matrix<samples::Scalar,N,N> Matrix;

    /// Return the number of dimensions
    int getDimension() const { return N; }

    /// Weight of this distribution in the mixture.
    Scalar weight;

    //@{
    /// Get/set the location parameter (mean/median/mode) of this component.
    Vector getMu() const { return _mu; }
    void setMu(Vector const & mu) { _mu = mu; }
    //@}

    //@{
    /**
     *  @brief Get/set the shape/size parameter
     *
     *  For the Gaussian distribution, this is simply the covariance matrix.
     *  For the Student's T distribution with df > 2, covariance = sigma * df / (df - 2);
     *  for df <= 2, the Student's T distribution has infinite variance, but is still a
     *  valid distribution.
     */
    Matrix getSigma() const { return _sigmaLLT.reconstructedMatrix(); }
    void setSigma(Matrix const & sigma) {
        _sigmaLLT.compute(sigma);
        _sqrtDet = _sigmaLLT.matrixLLT().diagonal().prod();
    }
    //@}

    /// Project the distribution onto the given dimension (marginalize over all others)
    MixtureComponent<1> project(int dim) const;

    /// Project the distribution onto the given dimensions (marginalize over all others)
    MixtureComponent<2> project(int dim1, int dim2) const;

    /// Default-construct a mixture component with weight=1, mu=0, sigma=identity.
    MixtureComponent() : weight(1.0), _sqrtDet(1.0), _mu(Vector::Zero()), _sigmaLLT(Matrix::Identity()) {}

    /// Default-construct a mixture component with the given parameters.
    MixtureComponent(Scalar weight_, Vector const & mu, Matrix const & sigma) :
        weight(weight_), _mu(mu), _sigmaLLT(sigma)
    {
        _sqrtDet = _sigmaLLT.matrixLLT().diagonal().prod();
    }

#ifndef SWIG
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    friend std::ostream & operator<< <>(std::ostream & os, MixtureComponent<N> const & self);
#endif

private:

    template <int M> friend class Mixture;

    void _stream(std::ostream & os, int offset=0) const;

    Scalar _sqrtDet;
    Vector _mu;
    Eigen::LLT<Matrix> _sigmaLLT;
};

/**
 *  @brief Helper class used to define restrictions to the form of the component parameters in
 *         Mixture::updateEM.
 *
 *  The base class implementation does not apply any restrictions.
 */
template <int N>
class MixtureUpdateRestriction {
public:
    typedef samples::Scalar Scalar;
    typedef typename MixtureComponent<N>::Vector Vector;
    typedef typename MixtureComponent<N>::Matrix Matrix;

    virtual void restrictMu(Vector & mu) const {}

    virtual void restrictSigma(Matrix & sigma) const {}

    virtual ~MixtureUpdateRestriction() {}
};

template <int N>
class Mixture :
    public afw::table::io::PersistableFacade< Mixture<N> >,
    public afw::table::io::Persistable
{
public:
    typedef samples::Scalar Scalar;
    typedef Eigen::Matrix<samples::Scalar,N,1> Vector;
    typedef Eigen::Matrix<samples::Scalar,N,N> Matrix;
    typedef MixtureComponent<N> Component;
    typedef MixtureUpdateRestriction<N> UpdateRestriction;
    typedef std::vector< Component, Eigen::aligned_allocator<Component> > ComponentList;
    typedef typename ComponentList::iterator iterator;
    typedef typename ComponentList::const_iterator const_iterator;

    //@{
    /**
     *  @brief Iterator and indexed access to components
     *
     *  While mutable iterators and accessors are provided, any modifications to
     *  the component weights should be followed by a call to normalize(), as
     *  other member functions will not work properly if the mixture is not
     *  normalized.
     */
    iterator begin() { return _components.begin(); }
    iterator end() { return _components.end(); }

    const_iterator begin() const { return _components.begin(); }
    const_iterator end() const { return _components.end(); }

    Component & operator[](std::size_t i) { return _components[i]; }
    Component const & operator[](std::size_t i) const { return _components[i]; }
    //@}

    /// Return the number of components
    std::size_t size() const { return _components.size(); }

    /// Project the distribution onto the given dimension (marginalize over all others)
    Mixture<1> project(int dim) const;

    /// Project the distribution onto the given dimensions (marginalize over all others)
    Mixture<2> project(int dim1, int dim2) const;

    /// Return the number of dimensions
    int getDimension() const { return N; }

    /// Iterate over all components, rescaling their weights so they sum to one.
    void normalize();

    /**
     *  @brief Iterate over all components, removing those with weight less than or equal to threshold
     *
     *  The weights will be normalized them if any are removed.
     *
     *  @return the number of components removed.
     */
    std::size_t clip(Scalar threshold=0.0);

    /// Get the number of degrees of freedom in the component Student's T distributions (inf=Gaussian)
    Scalar getDegreesOfFreedom() const { return _df; }

    /// Set the number of degrees of freedom in the component Student's T distributions (inf=Gaussian)
    void setDegreesOfFreedom(Scalar df=std::numeric_limits<Scalar>::infinity());

    /**
     *  @brief Evaluate the probability density at the given point for the given component distribution.
     *
     *  This evaluates the probability of a single component, including the current weight of that component.
     */
    template <typename Derived>
    Scalar evaluate(Component const & component, Eigen::MatrixBase<Derived> const & x) const {
        Scalar z = _computeZ(component, x);
        return component.weight * _evaluate(z) / component._sqrtDet;
    }

    /**
     *  @brief Evaluate the mixture distribution probability density function (PDF) at the given points
     *
     *  @param[in] x       point to evaluate, as an Eigen expression, shape=(N)
     */
    template <typename Derived>
    Scalar evaluate(Eigen::MatrixBase<Derived> const & x) const {
        Scalar p = 0.0;
        for (const_iterator i = begin(); i != end(); ++i) {
            p += evaluate(*i, x);
        }
        return p;
    }

    /**
     *  @brief Evaluate the mixture distribution probability density function (PDF) at the given points
     *
     *  @param[in] x       array of points, shape=(numSamples, N)
     *  @param[out] p      array of probability values, shape=(numSamples,)
     */
    void evaluate(ndarray::Array<Scalar const,2,1> const & x, ndarray::Array<Scalar,1,0> const & p) const;

    /**
     *  @brief Draw random variates from the mixture distribution.
     *
     *  @param[in,out] rng random number generator
     *  @param[out] x      array of points, shape=(numSamples, N)
     */
    void draw(afw::math::Random & rng, ndarray::Array<Scalar,2,1> const & x) const;

    /**
     *  @brief Perform an Expectation-Maximization step, updating the component parameters to match
     *         the given weighted samples.
     *
     *  @param[in] x       array of variables, shape=(numSamples, N)
     *  @param[in] w       array of weights, shape=(numSamples,)
     *  @param[in] restriction   Functor used to restrict the form of the updated mu and sigma
     */
    void updateEM(
        ndarray::Array<Scalar const,2,1> const & x,
        ndarray::Array<Scalar const,1,1> const & w,
        UpdateRestriction const & restriction=UpdateRestriction()
    );

    /**
     *  @brief Perform an Expectation-Maximization step, updating the component parameters to match
     *         the given unweighted samples.
     *
     *  @param[in] x       array of variables, shape=(numSamples, N)
     *  @param[in] restriction   Functor used to restrict the form of the updated mu and sigma
     */
    void updateEM(
        ndarray::Array<Scalar const,2,1> const & x,
        UpdateRestriction const & restriction=UpdateRestriction()
    );

    /**
     *  @brief Construct a mixture model.
     *
     *  @param[in] df            Number of degrees of freedom for component Student's T distributions
     *                           (inf=Gaussian)
     *  @param[in] components    List of components; will be emptied on return.
     *
     *  The components will be automatically normalized after construction.
     */
    explicit Mixture(ComponentList & components, Scalar df=std::numeric_limits<Scalar>::infinity());

#ifndef SWIG
    friend std::ostream & operator<< <>(std::ostream & os, Mixture<N> const & self);
#endif

    virtual bool isPersistable() const { return true; }

protected:

    virtual std::string getPythonModule() const { return "lsst.meas.multifit"; }

    virtual std::string getPersistenceName() const;

    virtual void write(OutputArchiveHandle & handle) const;

private:

    template <typename Derived>
    Scalar _computeZ(Component const & component, Eigen::MatrixBase<Derived> const & x) const {
        _workspace = x - component._mu;
        component._sigmaLLT.matrixL().solveInPlace(_workspace);
        return _workspace.squaredNorm();
    }

    Scalar _evaluate(Scalar z) const;

    void _stream(std::ostream & os) const;

    bool _isGaussian;
    Scalar _df;
    Scalar _norm;
    mutable Vector _workspace;
    ComponentList _components;
};

template <int N>
inline std::ostream & operator<<(std::ostream & os, MixtureComponent<N> const & self) {
    self._stream(os);
    return os;
}

template <int N>
inline std::ostream & operator<<(std::ostream & os, Mixture<N> const & self) {
    self._stream(os);
    return os;
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Mixture_h_INCLUDED
