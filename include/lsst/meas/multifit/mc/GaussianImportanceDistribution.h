// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_MC_GaussianImportanceDistribution
#define LSST_MEAS_MULTIFIT_MC_GaussianImportanceDistribution

#include "lsst/meas/multifit/mc/ImportanceDistribution.h"

namespace lsst { namespace meas { namespace multifit { namespace mc {

class GaussianImportanceDistribution : public ImportanceDistribution {
public:
    typedef boost::shared_ptr<GaussianImportanceDistribution> Ptr;

    /**
     *  @brief Draw a set of parameter vector from the distribution and evaluate the distribution
     *         at those points.
     *
     *  @param[in]   engine      Generic random number generator.
     *  @param[out]  parameters  (sample size)x(parameter count) array to fill with vectors drawn
     *                           from the distribution.
     *  @param[out]  importance  Density of the distribution.
     */
    virtual void draw(
        Random & engine,
        lsst::ndarray::Array<double,2,2> const & parameters,
        lsst::ndarray::Array<double,1,1> const & importance
    ) const;

    /**
     *  @brief Evaluate the distribution at the given parameter values.
     *
     *  @param[in]   parameters  (sample size)x(parameter count) array to evaluate at.
     *  @param[out]  output      (sample size) output array; the density times the given 
     *                           multiplicative factor is added to this array.
     *  @param[in]   factor      Multiplicative factor.
     */
    virtual void evaluate(
        lsst::ndarray::Array<double const,2,2> const & parameters,
        lsst::ndarray::Array<double,1,1> const & output,
        double factor = 1.0
    ) const;

    /**
     *  @brief Return a new distribution that has been modified match a set of weighted samples.
     *
     *  This will generally be used by adaptive importance sampling methods, and most
     *  operations will match moments or minimize the Kullback-Leibler divergence.
     */
    virtual ImportanceDistribution::Ptr adapt(
        lsst::ndarray::Array<double const,2,1> const & parameters,
        lsst::ndarray::Array<double const,1,1> const & weights
    ) const;

    static Ptr make(
        lsst::ndarray::Array<double const,1,1> const & mean,
        lsst::ndarray::Array<double const,2,2> const & sigma
    );

    static Ptr make(
        lsst::ndarray::Array<double const,1,1> const & mean,
        lsst::ndarray::Array<double const,1,1> const & sigmaDiagonal
    );

private:

    GaussianImportanceDistribution(Eigen::VectorXd & mean, Eigen::MatrixXd & sigma);

    double _normalization;
    Eigen::VectorXd _mean;
    Eigen::MatrixXd _factor;
    mutable Eigen::VectorXd _workspace;
};

}}}} // namespace lsst::meas::multifit::mc

#endif // !LSST_MEAS_MULTIFIT_MC_GaussianImportanceDistribution
