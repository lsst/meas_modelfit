#include "lsst/meas/multifit/ImportanceSampleTable.h"
#include "lsst/pex/exceptions.h"
#include "boost/format.hpp"

namespace lsst { namespace meas { namespace multifit {

ImportanceSampleTable ImportanceSampleTable::makeIterationTable(int n) const {
    return ImportanceSampleTable(*this, getIteration(n));
}

ImportanceSampleTable ImportanceSampleTable::makeLastIterationTable() const {
    return ImportanceSampleTable(*this, getLastIteration());
}

void ImportanceSampleTable::run(
    int size, Random & random, BaseDistribution const & importance,
    BaseEvaluator::Ptr const & evaluator
) {
    Evaluation evaluation(evaluator);
    Editor & editor = static_cast<Editor&>(_edit());
    editor.run(size, random, importance, evaluation);
}

void ImportanceSampleTable::run(
    int size, Random & random, BaseDistribution const & importance,
    BaseEvaluator::Ptr const & evaluator, BaseDistribution const & prior
) {
    Evaluation evaluation(evaluator); // TODO
    Editor & editor = static_cast<Editor&>(_edit());
    editor.run(size, random, importance, evaluation);
}

void ImportanceSampleTable::run(
    int size, Random & random, AlgorithmEnum algorithm,
    BaseEvaluator::Ptr const & evaluator
) {
    Evaluation evaluation(evaluator); // TODO
    Editor & editor = static_cast<Editor&>(_edit());
    editor.run(size, random, algorithm, evaluation);
}

void ImportanceSampleTable::run(
    int size, Random & random, AlgorithmEnum algorithm,
    BaseEvaluator::Ptr const & evaluator, BaseDistribution const & prior
) {
    Evaluation evaluation(evaluator);
    Editor & editor = static_cast<Editor&>(_edit());
    editor.run(size, random, algorithm, evaluation);
}

ImportanceSampleTable::ImportanceSampleTable(
    int capacity, int dimensionality, int nestedDimensionality, NestedMatrixType nestedMatrixType
) :
    NestedSampleTable(capacity, dimensionality, nestedDimensionality, nestedMatrixType),
    _iterations(),
    _objective(ndarray::allocate(capacity)),
    _importance(ndarray::allocate(capacity))  
{}

ImportanceSampleTable::ImportanceSampleTable(ImportanceSampleTable const & other) :
    NestedSampleTable(other),
    _iterations(other._iterations),
    _objective(other._objective),
    _importance(other._importance)
{}

ImportanceSampleTable::ImportanceSampleTable(
    ImportanceSampleTable const & other, Iteration const & iteration
) :
    NestedSampleTable(other, iteration.start, iteration.stop),
    _iterations(),
    _objective(other._objective[ndarray::view(iteration.start, iteration.stop)]),
    _importance(other._importance[ndarray::view(iteration.start, iteration.stop)])
{
    _iterations.push_back(Iteration(0, iteration.stop - iteration.start, iteration.distribution));
}

SampleTable::Ptr ImportanceSampleTable::_clone() const {
    return boost::make_shared<ImportanceSampleTable>(*this);
}

void ImportanceSampleTable::copyForEdit(int capacity) {
    NestedSampleTable::copyForEdit(capacity);
    copyArrayForEdit(_objective, capacity);
    copyArrayForEdit(_importance, capacity);
}

SampleTable::Editor::Ptr ImportanceSampleTable::makeEditor() {
    return boost::make_shared<Editor>(this);
}

void ImportanceSampleTable::Editor::run(
    int size, Random & random, BaseDistribution const & importance, Evaluation & evaluation
) {
    Iteration iteration(getSize(), size + getSize(), importance.clone());
    if (iteration.stop > getTable().getCapacity()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            (boost::format("New table size (%d) would be larger than capacity (%d).")
             % iteration.stop % getTable().getCapacity()).str()
        );
    }
    for (int n = iteration.start; n != iteration.stop; ++n) {
        iteration.distribution->draw(random, getParameters()[n]);
        evaluation.update(getParameters()[n]);
        getImportance()[n] = iteration.distribution->evaluate(evaluation.getParameters());
        // TODO
    }
    getSize() = iteration.stop;
}

void ImportanceSampleTable::Editor::run(
    int size, Random & random, AlgorithmEnum algorithm, Evaluation & evaluation
) {
    if (getTable().getIterationCount() == 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Adaptive importance sampling algorithms cannot be run unless "
            "at least one iteration already exists."
        );
    }
    Iteration const & last = getTable().getLastIteration();
    BaseDistribution::Ptr importance = last.distribution->clone();
    Iteration iteration(getSize(), size + getSize(), importance);
    if (iteration.stop > getTable().getCapacity()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            (boost::format("New table size (%d) would be larger than capacity (%d).")
             % iteration.stop % getTable().getCapacity()).str()
        );
    }
    importance->updateFromSamples(
        getParameters()[ndarray::view(last.start, last.stop)],
        getWeights()[ndarray::view(last.start, last.stop)]
    );
    for (int n = iteration.start; n != iteration.stop; ++n) {
        iteration.distribution->draw(random, getParameters()[n]);
        evaluation.update(getParameters()[n]);
        // TODO
    }
    getSize() = iteration.stop;
}

}}} // namespace lsst::meas::multifit
