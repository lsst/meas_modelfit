#include "lsst/meas/multifit/mc/SampleTable.h"
#include "lsst/pex/exceptions.h"
#include "boost/format.hpp"

namespace lsst { namespace meas { namespace multifit { namespace mc {

SampleTable::SampleTable(int capacity, int parameterCount) :
    _size(0),
    _parameters(ndarray::allocate(capacity, parameterCount)),
    _weights(ndarray::allocate(capacity))
{}

SampleTable::SampleTable(SampleTable const & other) :
    _size(other._size),
    _editor(other._editor),
    _parameters(other._parameters),
    _weights(other._weights)
{}

SampleTable::SampleTable(SampleTable const & other, int start, int stop) :
    _size(other._size),
    _editor(other._editor),
    _parameters(other._parameters[ndarray::view(start, stop)]),
    _weights(other._weights[ndarray::view(start, stop)])
{}

SampleTable::Editor & SampleTable::_edit(int capacity) {
    if (!_editor) {
        _editor = makeEditor();
    } else if (!_editor.unique()) {
        copyForEdit(capacity);
        _editor = makeEditor();
    } else if (capacity != getCapacity()) {
        copyForEdit(capacity);
    }
    return *_editor;
}

void SampleTable::copyForEdit(int capacity) {
    if (capacity < _size) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Requested capacity (%d) is smaller than current size (%d).")
             % capacity % _size).str()
        );
    }
    copyArrayForEdit(_parameters, capacity);
    copyArrayForEdit(_weights, capacity); 
}

}}}} // namespace lsst::meas::multifit::mc
