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
    int size, BaseEvaluator const & evaluator, Random & random, BaseDistribution const & distribution
) {
    Editor & editor = static_cast<Editor&>(_edit());
    editor.run(size, evaluator, random, distribution);
}

void ImportanceSampleTable::run(
    int size, BaseEvaluator const & evaluator, Random & random, AlgorithmEnum algorithm
) {
    Editor & editor = static_cast<Editor&>(_edit());
    editor.run(size, evaluator, random, algorithm);
}

ImportanceSampleTable::ImportanceSampleTable(
    int capacity, int dimensionality, int nestedDimensionality, NestedMatrixType nestedMatrixType
) :
    NestedSampleTable(capacity, dimensionality, nestedDimensionality, nestedMatrixType),
    _iterations(),
    _target(ndarray::allocate(capacity)),
    _importance(ndarray::allocate(capacity))  
{}

ImportanceSampleTable::ImportanceSampleTable(ImportanceSampleTable const & other) :
    NestedSampleTable(other),
    _iterations(other._iterations),
    _target(other._target),
    _importance(other._importance)
{}

ImportanceSampleTable::ImportanceSampleTable(
    ImportanceSampleTable const & other, Iteration const & iteration
) :
    NestedSampleTable(other, iteration.start, iteration.stop),
    _iterations(),
    _target(other._target[ndarray::view(iteration.start, iteration.stop)]),
    _importance(other._importance[ndarray::view(iteration.start, iteration.stop)])
{
    _iterations.push_back(Iteration(0, iteration.stop - iteration.start, iteration.distribution));
}

SampleTable::Ptr ImportanceSampleTable::_clone() const {
    return boost::make_shared<ImportanceSampleTable>(*this);
}

void ImportanceSampleTable::copyForEdit(int capacity) {
    NestedSampleTable::copyForEdit(capacity);
    copyArrayForEdit(_target, capacity);
    copyArrayForEdit(_importance, capacity);
}

SampleTable::Editor::Ptr ImportanceSampleTable::makeEditor() {
    return boost::make_shared<Editor>(this);
}

void ImportanceSampleTable::Editor::run(
    int size, BaseEvaluator const & evaluator, Random & random, BaseDistribution const & distribution
) {
    Iteration iteration(getSize(), size + getSize(), distribution.clone());
    if (iteration.stop > getTable().getCapacity()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            (boost::format("New table size (%d) would be larger than capacity (%d).")
             % iteration.stop % getTable().getCapacity()).str()
        );
    }
    for (int n = iteration.start; n != iteration.stop; ++n) {
        iteration.distribution->draw(random, getParameters()[n]);
        getWeights()[n] = iteration.distribution->evaluate(getParameters()[n]);
    }
    // TODO
}

}}} // namespace lsst::meas::multifit
