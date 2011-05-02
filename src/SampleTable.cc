#include "lsst/meas/multifit/SampleTable.h"
#include "lsst/pex/exceptions.h"
#include <boost/format.hpp>

namespace lsst { namespace meas { namespace multifit {

SampleTable::SampleTable(SampleTable const & other) :
    _size(other._size),
    _editor(other._editor),
    _parameters(other._parameters),
    _weights(other._weights)
{}

SampleTable::Editor & SampleTable::_edit() {
    if (!_editor) {
        _editor = makeEditor();
    } else if (!_editor.unique()) {
        copyForEdit(getCapacity());
        _editor = makeEditor();
    }
    return *_editor;
}

SampleTable::Editor & SampleTable::_reserve(int capacity) {
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
    lsst::ndarray::Array<double,2,2> newParameters(ndarray::allocate(capacity, getDimensionality()));
    lsst::ndarray::Array<double,1,1> newWeights(ndarray::allocate(capacity));
    newParameters[ndarray::view(0, _size)] = getParameters();
    newWeights[ndarray::view(0, _size)] = getWeights();
    _parameters = newParameters;
    _weights = newWeights;
}

void SampleTableEditor::_append(ndarray::Array<double const,1,1> const & parameters, double weight) {
    if (_table->getSize() >= _table->getCapacity()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Table is at capacity; cannot append additional records."
        );
    }
    _table->_parameters[_table->_size] = parameters;
    _table->_weights[_table->_size] = weight;
    ++_table->_size;
}

}}} // namespace lsst::meas::multifit
