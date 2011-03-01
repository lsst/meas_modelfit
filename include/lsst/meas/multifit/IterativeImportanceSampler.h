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

#ifndef LSST_MEAS_MULTIFIT_IterativeImportanceSampler
#define LSST_MEAS_MULTIFIT_IterativeImportanceSampler

#include "lsst/ndarray.h"
#include "lsst/ndarray/tables.h"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/random/uniform_real.hpp>

#include <vector>
#include <list>

namespace lsst { namespace meas { namespace multifit {

struct ImportanceSample {
    typedef boost::fusion::vector< 
        ndarray::tables::Field<double>,
        ndarray::tables::Field<double>,
        ndarray::tables::Field<double>,
        ndarray::tables::Field<double,1>,
        ndarray::tables::Field<double,1>,
        ndarray::tables::Field<double,2>
        > FieldSequence;

    static ndarray::tables::Index<0> const proposal;
    static ndarray::tables::Index<1> const target;
    static ndarray::tables::Index<2> const weight;
    static ndarray::tables::Index<3> const parameters;
    static ndarray::tables::Index<4> const coefficients;
    static ndarray::tables::Index<4> const fisher;

    typedef ndarray::tables::Table<Sample> Table;
    typedef ndarray::tables::Layout<Sample> Layout;
};

class MixtureDistribution {
public:

    typedef boost::mt19937 RandomEngine;

    class Component {
    public:

        Component(
            double normalization,
            Eigen::VectorXd const & mu,
            Eigen::MatrixXd const & sigma
        );

        double getNormalization() const { return _normalization; }

        Eigen::VectorXd const & getMu() const;

        Eigen::MatrixXd const & getSigma() const;

    private:
        double _normalization;
        Eigen::VectorXd _mu;
        Eigen::MatrixXd _sigma;
    };

    typedef std::list<Component> ComponentList;

    ComponentList const & getComponents() const { return _components; }

    ImportanceSample::Table draw(int size);

    MixtureDistribution(
        RandomEngine const & engine,
        ComponentList const & components
    );

private:
    typedef boost::variate_generator< RandomEngine &, boost::uniform_real<double> > UniformGenerator;
    typedef boost::variate_generator< RandomEngine &, boost::normal_distribution<double> > NormalGenerator;

    RandomEngine _randomEngine;
    UniformGenerator _randomUniform;
    NormalGenerator _randomNormal;
    ComponentList _components;
};


class IterativeImportanceSampler {
public:

    int getIterationCount() const { return _samples.size(); }

    Table const & getTable(int n) const { return _samples[n]; }

    BaseEvaluator::Ptr getEvaluator() const;

    void run(int size);

private:
    BaseEvaluator::Ptr _evaluator;
    MixtureDistribution _proposal;
    std::vector<Table> _samples;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_IterativeImportanceSampler
