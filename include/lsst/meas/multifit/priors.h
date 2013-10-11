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
     *     P(D|\alpha,\theta) = A e^{-g^T\alpha - \frac{1}{2}\alpha^T F \alpha}
     *  @f]
     *  The normalization @f$A@f$ can be brought outside the integral as a constant to be added
     *  to the return value, so it is not passed as an argument to this function.
     *
     *  @param[in]  gradient     Gradient of the -log likelihood in @f$\alpha@f$ at fixed @f$\theta@f$;
     *                           the vector @f$g@f$ in the equation above.
     *  @param[in]  fisher       Second derivatives of of the -log likelihood in @f$\alpha@f$ at fixed
     *                           @f$\theta@f$; the matrix @f$F@f$ in the equation above.
     */
    virtual Scalar marginalize(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear
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

    /// @copydoc Prior::evaluate
    virtual Scalar marginalize(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear
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
