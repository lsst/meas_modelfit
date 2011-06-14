#include "lsst/meas/multifit/mc/NestedSampleTable.h"

namespace lsst { namespace meas { namespace multifit { namespace mc {

NestedSampleTable::NestedSampleTable(
    int capacity, int nestedSize, int parameterCount, int coefficientCount
) :
    SampleTable(capacity, parameterCount),
    _nestedWeights(ndarray::allocate(capacity, nestedSize)),
    _coefficients(ndarray::allocate(capacity, nestedSize, coefficientCount))
{}

NestedSampleTable::NestedSampleTable(NestedSampleTable const & other) :
    SampleTable(other),
    _nestedWeights(other._nestedWeights),
    _coefficients(other._coefficients)
{}

NestedSampleTable::NestedSampleTable(NestedSampleTable const & other, int start, int stop) :
    SampleTable(other, start, stop),
    _nestedWeights(other._nestedWeights[ndarray::view(start, stop)]),
    _coefficients(other._coefficients[ndarray::view(start, stop)])
{}

void NestedSampleTable::copyForEdit(int capacity) {
    SampleTable::copyForEdit(capacity);
    copyArrayForEdit(_nestedWeights, capacity);
    copyArrayForEdit(_coefficients, capacity); 
}

}}}} // namespace lsst::meas::multifit::mc
