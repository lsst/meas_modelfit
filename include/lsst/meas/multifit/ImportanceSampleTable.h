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

#ifndef LSST_MEAS_MULTIFIT_ImportanceSampleTable
#define LSST_MEAS_MULTIFIT_ImportanceSampleTable

#include "lsst/meas/multifit/NestedSampleTable.h"
#include "lsst/meas/multifit/BaseDistribution.h"
#include "lsst/meas/multifit/BaseEvaluator.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A SampleTable for use in adaptive importance sampling algorithms.
 */
class ImportanceSampleTable : public NestedSampleTable {
public:

    typedef boost::shared_ptr<ImportanceSampleTable> Ptr;
    typedef boost::shared_ptr<ImportanceSampleTable const> ConstPtr;

    /**
     *  @brief Enum to select the algorithm used in iteratively updating distributions.
     */
    enum AlgorithmEnum {
        AIS, /**< Adaptive Importance Sampling; use the last group of samples to update the distribution.
              *   See Cappe et al 2008, arXiv:0710.4242v4
              */
        AMIS /**< Adaptive Multiple Importance Sampling; use all samples to update the distribution.
              *   See Cornuet et al 2010, arXiv:0907.1254v3
              */
    };

#ifndef SWIG

    /**
     *  @brief A simple struct to represent one iteration of an adaptive importance sampling run
     *         within a full ImportanceSampleTable.
     *
     *  Iterations may overlap if the AMIS algorithm is used; an AMIS run will create an iteration
     *  that contains the entire previous iteration.
     */
    struct Iteration {
        int start;
        int stop;
        BaseDistribution::ConstPtr distribution;

        Iteration(int start_, int stop_, BaseDistribution::ConstPtr const & distribution_) :
            start(start_), stop(stop_), distribution(distribution_)
        {}
    };

#endif

    /**
     *  @brief The unnormalized value of the target distribution at the sample points.
     */
    lsst::ndarray::Array<double const,1,1> getTarget() const {
        return _target[ndarray::view(0, getSize())];
    }

    /**
     *  @brief The density of the distribution the sample was drawn from, multiplied by
     *         the total number of samples.
     */
    lsst::ndarray::Array<double const,1,1> getImportance() const {
        return _importance[ndarray::view(0, getSize())];
    }

    /// @brief Return the number of iterations.
    int getIterationCount() const { return _iterations.size(); }

#ifndef SWIG
    /// @brief Return the indices and importance distribution for the given iteration.
    Iteration const & getIteration(int n) const { return _iterations[n]; }

    /// @brief Return the last iteration.
    Iteration const & getLastIteration() const { return _iterations.back(); }
#endif

    /// @brief Get the subset of the table corresponding to a particular iteration as a new table.
    ImportanceSampleTable makeIterationTable(int n) const;

    /// @brief Get the subset of the table corresponding to the last iteration as a new table.
    ImportanceSampleTable makeLastIterationTable() const;

    /// @brief Copy the table (shallow with copy-on-write).
    Ptr clone() const { return boost::static_pointer_cast<ImportanceSampleTable>(_clone()); }
  
    /**
     *  @brief Set the capacity of the table by reallocating all arrays and return an Editor.
     *
     *  The capacity may not be less than the current size.
     */
    void reserve(int capacity) { _edit(capacity); }

    /**
     *  @brief Add an iteration of standard importance sampling using the given distribution.
     */
    void run(int size, BaseEvaluator const & evaluator, Random & random, 
             BaseDistribution const & distribution);
    
    /**
     *  @brief Run an adaptive importance sampling algorithm, using an updated version
     *         of the last iteration's distribution.
     *
     *  The table must have at least one iteration before this overload can be called.
     *
     *  Note that both algorithms will use only the last Iteration's samples to update
     *  the distribution, but the AMIS algorithm will include those samples in its own
     *  Iteration, so a sequence of AMIS runs will use all samples.  This allows
     *  different algorithms to be used for different iterations.
     *
     *  While the user can use a different evaluator for different iterations, they should
     *  generate similar likelihoods (up to the normization) for the procedure to be robust.
     *  The AIS algorithm must be used when the evaluator is changed.
     *
     *  In all cases, only the last iteration can be considered a fair sample.
     */
    void run(int size, BaseEvaluator const & evaluator, Random & random, AlgorithmEnum algorithm);

    /// @brief Construct with zero size and finite capacity.
    explicit ImportanceSampleTable(
        int capacity, int dimensionality, int nestedDimensionality, NestedMatrixType nestedMatrixType
    );

    /// @brief Copy constructor.
    ImportanceSampleTable(ImportanceSampleTable const & other);

protected:

#ifndef SWIG
    class Editor : public NestedSampleTable::Editor {
    public:
        void run(int size, BaseEvaluator const & evaluator, Random & random, 
                 BaseDistribution const & distribution);
        void run(int size, BaseEvaluator const & evaluator, Random & random, AlgorithmEnum algorithm);

        explicit Editor(ImportanceSampleTable * table) : NestedSampleTable::Editor(table) {}

    protected:

        lsst::ndarray::Array<double,1,1> getTarget() {
            return getTable()._target;
        }

        lsst::ndarray::Array<double,1,1> getImportance() {
            return getTable()._importance;
        }

        ImportanceSampleTable & getTable() {
            return static_cast<ImportanceSampleTable&>(NestedSampleTable::Editor::getTable());
        }
    };
#endif

    ImportanceSampleTable(ImportanceSampleTable const & other, Iteration const & iteration);

    virtual SampleTable::Ptr _clone() const;
    
    virtual void copyForEdit(int capacity);

    virtual SampleTable::Editor::Ptr makeEditor();

private:
    std::vector<Iteration> _iterations;
    lsst::ndarray::Array<double,1,1> _target;
    lsst::ndarray::Array<double,1,1> _importance;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_ImportanceSampleTable
