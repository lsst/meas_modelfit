#include "lsst/meas/multifit/mc/ImportanceSampleTable.h"
#include "lsst/pex/exceptions.h"
#include "boost/format.hpp"
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace mc {

void ImportanceSampleTable::run(int size, Random & random) {
    Editor & editor = static_cast<Editor&>(_edit(std::max(size + getSize(), getCapacity())));
    editor.run(size, random);
}
ImportanceSampleTable::ImportanceSampleTable(
    ImportanceDistribution::Ptr const & distribution,
    BaseEvaluator::Ptr const & evaluator,
    int nestedSize, int capacity
) :
    NestedSampleTable(capacity, nestedSize, evaluator->getParameterCount(), evaluator->getCoefficientCount()),
    _iterations(),
    _distribution(distribution),
    _evaluator(evaluator),
    _objective(ndarray::allocate(capacity)),
    _importance(ndarray::allocate(capacity))
{}

ImportanceSampleTable::ImportanceSampleTable(ImportanceSampleTable const & other) :
    NestedSampleTable(other),
    _iterations(other._iterations),
    _distribution(other._distribution),
    _evaluator(other._evaluator),
    _objective(other._objective),
    _importance(other._importance)
{}

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
    int size, Random & random
) {
    int start = getSize();
    int stop = start + size;
    assert(stop <= getTable().getCapacity());
    if (!getTable()._iterations.empty()) {
        assert(start != 0);
        getTable()._distribution = getTable()._distribution->adapt(
            getParameters()[ndarray::view(0, start)],
            getWeights()[ndarray::view(0, start)]
        );
    }
    getTable()._distribution->draw(
        random, 
        getParameters()[ndarray::view(start, stop)],
        getImportance()[ndarray::view(start, stop)]
    );
    ndarray::EigenView<double,1,1> importance(getImportance()[ndarray::view(0, stop)]);
    if (!getTable()._iterations.empty()) {
        importance.segment(start, size) *= size;
        for (
            IterationVector::iterator i = getTable()._iterations.begin(); 
            i != getTable()._iterations.end();
            ++i
        ) {
            i->first->evaluate(
                getParameters()[ndarray::view(start, stop)],
                getImportance()[ndarray::view(start, stop)],
                i->second
            );
        }
        importance.segment(0, start) *= getSize();
        getTable()._distribution->evaluate(
            getParameters()[ndarray::view(0, start)],
            getImportance()[ndarray::view(0, start)],
            size
        );
        importance /= stop;
    }
    for (int n = start; n != stop; ++n) {
        getObjective()[n] = getTable()._evaluator->integrate(
            random,
            getCoefficients()[n],
            getNestedWeights()[n],
            getParameters()[n]
        );
    }
    ndarray::viewAsEigen(getWeights()).segment(0, stop) 
        = ndarray::viewAsEigen(getObjective()).segment(0, stop).cwise() / importance;
    getTable()._iterations.push_back(std::make_pair(getTable()._distribution, size));
    getSize() = stop;
}

}}}} // namespace lsst::meas::multifit::mc
