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

#ifndef LSST_MEAS_MULTIFIT_priors_h_INCLUDED
#define LSST_MEAS_MULTIFIT_priors_h_INCLUDED

#include "lsst/base.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/afw/math/Random.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/Mixture.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for Bayesian priors
 */
class Prior :
    public afw::table::io::PersistableFacade<Prior>,
    public afw::table::io::Persistable,
    private boost::noncopyable
{
public:

    std::string const & getTag() const { return _tag; }

    /**
     *  @brief Evaluate the prior at the given point in nonlinear and amplitude space.
     *
     *  @param[in]   nonlinear        Vector of nonlinear parameters
     *  @param[in]   amplitudes       Vector of linear parameters
     */
    virtual Scalar evaluate(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes
    ) const = 0;

    /**
     *  @brief Evaluate the derivatives of the prior at the given point in nonlinear and amplitude space.
     *
     *  Note that while the model is linear in the amplitudes, the prior is not necessarily
     *  linear in the amplitudes, so we do care about second derivatives w.r.t. amplitudes.
     *
     *  @param[in]   nonlinear           Vector of nonlinear parameters
     *  @param[in]   amplitudes          Vector of linear parameters
     *  @param[in]   nonlinearGradient   First derivative w.r.t. nonlinear parameters
     *  @param[in]   amplitudeGradient   First derivative w.r.t. linear parameters parameters
     *  @param[in]   nonlinearHessian    Second derivative w.r.t. nonlinear parameters
     *  @param[in]   amplitudeHessian    Second derivative w.r.t. linear parameters parameters
     *  @param[in]   crossHessian        Second derivative cross term of d(nonlinear)d(amplitudes);
     *                                   shape is [nonlinearDim, amplitudeDim].
     */
    virtual void evaluateDerivatives(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & nonlinearGradient,
        ndarray::Array<Scalar,1,1> const & amplitudeGradient,
        ndarray::Array<Scalar,2,1> const & nonlinearHessian,
        ndarray::Array<Scalar,2,1> const & amplitudeHessian,
        ndarray::Array<Scalar,2,1> const & crossHessian
    ) const = 0;

    /**
     *  @brief Return the -log amplitude integral of the prior*likelihood product.
     *
     *  If @f$\alpha@f$ are the amplitudes, @f$\theta@f$ are the nonlinear parameters, and @f$D@f$ is
     *  the data, then this method represents @f$P(\alpha,\theta)@f$ by computing
     *  @f[
     *    -\ln\left[\int\!P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha\right]
     *  @f]
     *  at fixed @f$\theta@f$.  Because @f$\alpha@f$ are linear parameters, @f$P(D|\alpha,\theta)@f$
     *  is Gaussian in @f$\alpha@f$, and because @f$\theta@f$ is fixed, it's usually convenient to
     *  think of the integral as:
     *  @f[
     *    -ln\left[P(\theta)\int\!P(D|\alpha,\theta)\,P(\alpha|\theta)\,d\alpha\right]
     *  @f]
     *  Thus, we marginalize the likelihood in @f$\alpha@f$ at fixed @f$\theta@f$, and then multiply
     *  by the prior on @f$\theta@f$.
     *
     *  We also assume the likelihood @f$P(D|\alpha,\theta)@f$ is Gaussian in @f$\alpha@f$, which is
     *  generally true because @f$\alpha@f$ defined such that the model is linear in them, and the
     *  noise on the data is generally Gaussian.  In detail, we represent the likelihood at fixed
     *  @f$\theta@f$ as
     *  @f[
     *     P(D|\alpha,\theta) = A e^{-g^T\alpha - \frac{1}{2}\alpha^T H \alpha}
     *  @f]
     *  The normalization @f$A@f$ can be brought outside the integral as a constant to be added
     *  to the return value, so it is not passed as an argument to this function.
     *
     *  @param[in]  gradient     Gradient of the -log likelihood in @f$\alpha@f$ at fixed @f$\theta@f$;
     *                           the vector @f$g@f$ in the equation above.
     *  @param[in]  hessian      Second derivatives of of the -log likelihood in @f$\alpha@f$ at fixed
     *                           @f$\theta@f$; the matrix @f$H@f$ in the equation above.
     *  @param[in]  nonlinear    The nonlinear parameters @f$\theta@f$.
     */
    virtual Scalar marginalize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear
    ) const = 0;

    /**
     *  @brief Draw a set of Monte Carlo amplitude vectors
     *
     *  This provides a Monte Carlo approach to extracting the conditional amplitude distribution that
     *  is integrated by the marginalize() method.
     *
     *  @param[in]  gradient     Gradient of the -log likelihood in @f$\alpha@f$ at fixed @f$\theta@f$.
     *  @param[in]  fisher       Second derivatives of of the -log likelihood in @f$\alpha@f$ at fixed
     *                           @f$\theta@f$.
     *  @param[in]  nonlinear    The nonlinear parameters @f$\theta@f$ at which we are evaluating
     *                           the conditional distribution @f$P(\alpha|\theta)@f$.
     *  @param[in,out]  rng      Random number generator.
     *  @param[out] amplitudes   The Monte Carlo sample of amplitude parameters @f$\alpha@f$.  The
     *                           number of rows sets the number of samples, while the number of
     *                           columns must match the dimensionality of @f$\alpha@f$.
     *  @param[out] weights      The weights of the Monte Carlo samples; should asymptotically average
     *                           to one.
     *  @param[in]  multiplyWeights  If true, multiply weight vector instead of overwriting it.
     */
    virtual void drawAmplitudes(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        afw::math::Random & rng,
        ndarray::Array<Scalar,2,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & weights,
        bool multiplyWeights=false
    ) const = 0;

    virtual ~Prior() {}

protected:

    explicit Prior(std::string const & tag="") : _tag(tag) {}

    virtual std::string getPythonModule() const { return "lsst.meas.multifit"; }

private:
    std::string _tag;
};

/**
 *  @brief A prior that's flat in amplitude parameters, and uses a Mixture for nonlinear parameters.
 */
class MixturePrior :
    public afw::table::io::PersistableFacade<MixturePrior>,
    public Prior
{
public:

    explicit MixturePrior(PTR(Mixture const) mixture, std::string const & tag="");

    /// @copydoc Prior::evaluate
    virtual Scalar evaluate(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes
    ) const;

    /// @copydoc Prior::evaluateDerivatives
    virtual void evaluateDerivatives(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & nonlinearGradient,
        ndarray::Array<Scalar,1,1> const & amplitudeGradient,
        ndarray::Array<Scalar,2,1> const & nonlinearHessian,
        ndarray::Array<Scalar,2,1> const & amplitudeHessian,
        ndarray::Array<Scalar,2,1> const & crossHessian
    ) const;

    /// @copydoc Prior::evaluate
    virtual Scalar marginalize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear
    ) const;

    /// @copydoc Prior::drawAmplitudes
    virtual void drawAmplitudes(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        afw::math::Random & rng,
        ndarray::Array<Scalar,2,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & weights,
        bool multiplyWeights=false
    ) const;

    /**
     *  @brief Return a MixtureUpdateRestriction appropriate for (e1,e2,r) data.
     *
     *  This restriction object can be used with Mixture<3>::updateEM() to create
     *  a mixture with a strictly isotropic ellipticity distribution.
     */
    static MixtureUpdateRestriction const & getUpdateRestriction();

    PTR(Mixture const) getMixture() const { return _mixture; }

    virtual bool isPersistable() const { return true; }

protected:

    virtual std::string getPersistenceName() const;

    virtual void write(OutputArchiveHandle & handle) const;

private:
    PTR(Mixture const) _mixture;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_priors_h_INCLUDED
