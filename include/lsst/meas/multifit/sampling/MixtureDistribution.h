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

#include <lsst/ndarray.h>
#include <Eigen/Core>

#include <boost/random/mersenne_twister.hpp>

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

    Eigen::VectorXd const & getMu() const { return _mu; }

    Eigen::Part<Eigen::MatrixXd,Eigen::LowerTriangular> const getSigma() const {
        return _sigma.part<Eigen::LowerTriangular>();
    }

private:

    friend class MixtureDistribution;

    double _normalization;
    Eigen::VectorXd _mu;
    Eigen::MatrixXd _sigma;
};

class MixtureDistribution {
public:

    typedef boost::mt19937 RandomEngine;
    typedef MixtureComponent Component;
    typedef std::vector<Component> ComponentList;

    explicit MixtureDistribution(ComponentList const & components);

    ComponentList const & getComponents() const { return _components; }

    void draw(
        lsst::ndarray::Array<double,2,2> const & points,
        RandomEngine & engine
    ) const;

    void evaluate(
        lsst::ndarray::Array<double,1,1> const & probability,
        lsst::ndarray::Array<double const,2,2> const & points
    ) const;

private:
    ComponentList _components;
};

}}}} // namespace lsst::meas::multifit::sampling

#endif // !LSST_MEAS_MULTIFIT_SAMPLING_MixtureDistribution
