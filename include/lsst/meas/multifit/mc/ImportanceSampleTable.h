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

#ifndef LSST_MEAS_MULTIFIT_MC_ImportanceSampleTable
#define LSST_MEAS_MULTIFIT_MC_ImportanceSampleTable

#include "lsst/meas/multifit/mc/NestedSampleTable.h"
#include "lsst/meas/multifit/BaseEvaluator.h"

namespace lsst { namespace meas { namespace multifit { namespace mc {

class ImportanceDistribution {
public:
    typedef boost::shared_ptr<ImportanceDistribution> Ptr;

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
    ) const = 0;

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
    ) const = 0;

    /**
     *  @brief Return a new distribution that has been modified match a set of importance or MCMC samples.
     *
     *  This will generally be used by adaptive importance sampling methods, and most
     *  operations will match moments or minimize the Kullback-Leibler divergence.
     */
    virtual Ptr adapt(
        lsst::ndarray::Array<double const,2,1> const & parameters,
        lsst::ndarray::Array<double const,1,1> const & weights
    ) const = 0;

};

/**
 *  @brief A SampleTable for use in adaptive importance sampling algorithms.
 */
class ImportanceSampleTable : public NestedSampleTable {
public:

    typedef boost::shared_ptr<ImportanceSampleTable> Ptr;
    typedef boost::shared_ptr<ImportanceSampleTable const> ConstPtr;

    /**
     *  @brief Fill the table using importance sampling.
     *
     *  @param[in] size        Number of samples to add to the table.
     *  @param[in] random      Random number generator.
     *
     *  If the current size of the table is nonzero, the "Adaptive Multiple Importance Sampling"
     *  algorithm of Cornuet et al 2010 (arXiv:0907.1254v3) is used to append additional samples
     *  and update the importance and weights values of existing samples.  The distribution
     *  will be updated to match the current samples before new samples are drawn.
     *
     *  If the final size of the table is less than the capacity, the array will be reallocated.
     *  Unlike std::vector, the capacity will not be increased beyond the new size automatically,
     *  so it is important to use reserve() or construct the table with the desired capacity in
     *  advance whenever possible.
     *
     *  To run non-cumulative adaptive importance sampling (c.f. Cappe et al 2008, arXiv:0710.4242v4)
     *  alternate calls to run() and reset(true).
     */
    void run(int size, Random & random);

    /**
     *  @brief Clear all records.
     *
     *  The table will only be reallocated and the capacity changed if it shares data with 
     *  another table (i.e. copy-on-write).
     *
     *  The evaluator may also be changed during the reset; passing in an empty evaluator
     *  causes the table to keep the current evaluator.  The new evaluator must have
     *  the same parameter and coefficient dimensions as the original one.
     */
    void reset(bool adapt=true, BaseEvaluator::Ptr const & evaluator=BaseEvaluator::Ptr());

    /**
     *  @brief Clear all records and change the importance distribution.
     *
     *  The table will only be reallocated and the capacity changed if it shares data with 
     *  another table (i.e. copy-on-write).
     *
     *  The evaluator may also be changed during the reset; passing in an empty evaluator
     *  causes the table to keep the current evaluator.  The new evaluator must have
     *  the same parameter and coefficient dimensions as the original one.
     */
    void reset(
        ImportanceDistribution::Ptr const & distribution,
        BaseEvaluator::Ptr const & evaluator=BaseEvaluator::Ptr()
    );

    /**
     *  @brief Set the capacity of the table by reallocating all arrays.
     *
     *  The capacity may not be less than the current size.
     */
    void reserve(int capacity) { _edit(capacity); }

    /**
     *  @brief The value of the objective function at each sample point.
     */
    lsst::ndarray::Array<double const,1,1> getObjective() const {
        return _objective[ndarray::view(0, getSize())];
    }

    /**
     *  @brief The density of the distribution the sample was drawn from.
     */
    lsst::ndarray::Array<double const,1,1> getImportance() const {
        return _importance[ndarray::view(0, getSize())];
    }

    /// @brief Copy the table (shallow with copy-on-write).
    Ptr clone() const { return boost::static_pointer_cast<ImportanceSampleTable>(_clone()); }

    /**
     *  @brief Return the distribution samples are drawn from.
     *
     *  If the table has been run cumulatively, this only returns the last distribution,
     *  not the pseudo-mixture distribution.
     */
    ImportanceDistribution::Ptr getDistribution() const { return _distribution; }

    /**
     *  @brief Set the distribution samples are drawn from.
     */
    void setDistribution(ImportanceDistribution::Ptr const & distribution) { _distribution = distribution; }

    /// @brief Return the evaluator that defines the objective function.
    BaseEvaluator::Ptr getEvaluator() const { return _evaluator; }

    /**
     *  @brief Construct a table with the given importance distribution, evaluator, and nested size.
     *
     *  @param[in]  distribution   Distribution to draw samples from.
     *  @param[in]  evaluator      Evaluator that defines the objective function.
     *  @param[in]  nestedSize     Number of coefficient samples to draw for each parameter sample.
     *  @param[in]  capacity       Anticipated number of (parameter) samples in the table.
     */
    ImportanceSampleTable(
        ImportanceDistribution::Ptr const & distribution,
        BaseEvaluator::Ptr const & evaluator,
        int nestedSize, int capacity=0
    );

    /// @brief Copy constructor.
    ImportanceSampleTable(ImportanceSampleTable const & other);

protected:

#ifndef SWIG
    class Editor : public NestedSampleTable::Editor {
    public:

        explicit Editor(ImportanceSampleTable * table) : NestedSampleTable::Editor(table) {}

        void run(int size, Random & random);

    protected:

        ndarray::Array<double,1,1> const & getObjective() {
            return getTable()._objective;
        }

        ndarray::Array<double,1,1> const & getImportance() {
            return getTable()._importance;
        }

        ImportanceSampleTable & getTable() {
            return static_cast<ImportanceSampleTable&>(NestedSampleTable::Editor::getTable());
        }
    };
#endif

    virtual SampleTable::Ptr _clone() const;
    
    virtual void copyForEdit(int capacity);

    virtual SampleTable::Editor::Ptr makeEditor();

private:
    typedef std::pair<ImportanceDistribution::Ptr,int> Iteration;
    typedef std::vector<Iteration> IterationVector;

    IterationVector _iterations;
    ImportanceDistribution::Ptr _distribution;
    BaseEvaluator::Ptr _evaluator;
    ndarray::Array<double,1,1> _objective;
    ndarray::Array<double,1,1> _importance;
};

}}}} // namespace lsst::meas::multifit::mc

#endif // !LSST_MEAS_MULTIFIT_MC_ImportanceSampleTable
