#include "lsst/meas/multifit/NestedSampleTable.h"

namespace lsst { namespace meas { namespace multifit {

NestedSampleTable::NestedSampleTable(
    int capacity, int dimensionality, int nestedDimensionality,
    NestedMatrixType nestedMatrixType
) :
    SampleTable(capacity, dimensionality),
    _nestedMatrixType(nestedMatrixType),
    _nestedAmplitudes(ndarray::allocate(capacity)),
    _nestedMeans(ndarray::allocate(capacity, nestedDimensionality)),
    _nestedMatrices(ndarray::allocate(capacity, nestedDimensionality, nestedDimensionality))
{}

NestedSampleTable::NestedSampleTable(NestedSampleTable const & other) :
    SampleTable(other),
    _nestedMatrixType(other._nestedMatrixType),
    _nestedAmplitudes(other._nestedAmplitudes),
    _nestedMeans(other._nestedMeans),
    _nestedMatrices(other._nestedMatrices)
{}

NestedSampleTable::NestedSampleTable(NestedSampleTable const & other, int start, int stop) :
    SampleTable(other, start, stop),
    _nestedMatrixType(other._nestedMatrixType),
    _nestedAmplitudes(other._nestedAmplitudes[ndarray::view(start, stop)]),
    _nestedMeans(other._nestedMeans[ndarray::view(start, stop)]),
    _nestedMatrices(other._nestedMatrices[ndarray::view(start, stop)])
{}

void NestedSampleTable::copyForEdit(int capacity) {
    SampleTable::copyForEdit(capacity);
    copyArrayForEdit(_nestedAmplitudes, capacity);
    copyArrayForEdit(_nestedMeans, capacity); 
    copyArrayForEdit(_nestedMatrices, capacity); 
}

}}} // namespace lsst::meas::multifit
