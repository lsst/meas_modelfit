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

#ifndef LSST_MEAS_MULTIFIT_SAMPLING_MixtureDistribution
#define LSST_MEAS_MULTIFIT_SAMPLING_MixtureDistribution

#include "lsst/ndarray.h"
#include "lsst/meas/multifit/sampling/Table.h"
#include "lsst/meas/multifit/sampling/RandomEngine.h"

#include <Eigen/Core>

#include <vector>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

class MixtureComponent {
public:

    MixtureComponent(
        double normalization,
        Eigen::VectorXd const & mu,
        Eigen::MatrixXd const & sigma
    ) : _normalization(normalization), _mu(mu), _sigma(sigma) {}

    double getNormalization() const { return _normalization; }

    /// @brief Mean vector.
    Eigen::VectorXd const & getMu() const { return _mu; }

    /// @brief Lower triangular Cholesky factor of the covariance matrix.
    Eigen::MatrixXd const & getSigma() const { return _sigma; }

private:

    friend class MixtureDistribution;

    double _normalization;
    Eigen::VectorXd _mu;
    Eigen::MatrixXd _sigma;
};

class MixtureDistribution {
public:

    typedef MixtureComponent Component;
    typedef std::vector<Component> ComponentList;

    explicit MixtureDistribution(ComponentList const & components, int dof=-1);

    ComponentList const & getComponents() const { return _components; }

    int getParameterSize() const { return _components.front().getMu().size(); }

    void draw(Table const & table, RandomEngine & engine) const;

    void update(ConstTable const & table);

    /**
     *  @brief Create a mixture distribution from a single parameter vector and covariance matrix.
     *
     *  The components of the returned mixture will be drawn a Gaussian defined by the given
     *  mean and the covariance scaled by the given fraction.
     *
     *  @param[in[ nComponents      Number of components in the mixture.
     *  @param[in] fraction         Fraction of covariance to scatter mean values with.
     *  @param[in] mean             Mean vector.
     *  @param[in] covariance       Covariance matrix.
     */
    static MixtureDistribution scatter(
        RandomEngine & engine,
        int nComponents, double fraction, int dof,
        lsst::ndarray::Array<double const,1,1> const & mean,
        lsst::ndarray::Array<double const,2,2> const & covariance
    );

private:
    int _dof;
    double _constant;
    ComponentList _components;
};

}}}} // namespace lsst::meas::multifit::sampling

#endif // !LSST_MEAS_MULTIFIT_SAMPLING_MixtureDistribution
