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

#ifndef LSST_MEAS_MULTIFIT_SampleTable
#define LSST_MEAS_MULTIFIT_SampleTable

#include "lsst/ndarray.h"
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace lsst { namespace meas { namespace multifit {

class SampleTableEditor;

/**
 *  A base class for tables of weighted samples of points in parameter space.
 *
 *  While one cannot append to the base SampleTable class, leaf classes will generally provide
 *  an append method that adds additional records.  The whole table is allocated in advance, however,
 *  and (unlike std::vector) one cannot append additional rows when the table is at capacity
 *  unless reserve() is called explicitly.
 */
class SampleTable {
public:

    typedef boost::shared_ptr<SampleTable> Ptr;
    typedef boost::shared_ptr<SampleTable const> ConstPtr;

    typedef SampleTableEditor Editor;

    /**
     *  @brief The immutable (samples)x(parameters) matrix of parameter values.
     */
    lsst::ndarray::Array<double const,2,2> getParameters() const {
        return _parameters[ndarray::view(0, _size)];
    }

    /**
     *  @brief The weights corresponding to each sample.
     */
    lsst::ndarray::Array<double const,1,1> getWeights() const {
        return _weights[ndarray::view(0, _size)];
    }

    /// @brief Return the number of parameters in the sample table.
    int getDimensionality() const { return _parameters.getSize<1>(); }

    /// @brief Return the number of samples in the table.
    int getSize() const { return _size; }

    /// @brief Return the number of samples that the table has been allocated to support.
    int getCapacity() const { return _parameters.getSize<0>(); }

    /**
     *  @brief Copy the table (shallow with copy-on-write).
     *
     *  Subclasses should reimplement, calling _clone() and casting the returned pointer
     *  to the appropriate subclass of SampleTable.
     */
    Ptr clone() const { return _clone(); }

    /**
     *  @brief Set the capacity of the table by reallocating all arrays and return an Editor.
     *
     *  The capacity may not be less than the current size.
     *
     *  Subclasses should reimplement, calling _reserve() and casting the returned reference
     *  to the appropriate subclass of Editor.
     */
    Editor & reserve(int capacity) { return _reserve(capacity); }

    /**
     *  @brief Return an interface class that allows changes to be made to the table.
     *
     *  Subclasses should reimplement, calling _edit() and casting the returned reference
     *  to the appropriate subclass of Editor.
     */
    Editor & edit() { return _edit(); }

    virtual ~SampleTable() {}

protected:

    friend class SampleTableEditor;

    SampleTable(SampleTable const & other);

    /**
     *  @Brief Implementation for edit(); moved here so subclasses can cast the public return value.
     *
     *  The copy-on-write mechanism works by checking the _editor data member.  If it is empty (as
     *  it will be for a newly-constructed, non-copied table) makeEditor() is called to construct
     *  a new Editor.  If it is unique, the existing Editor is returned.  If it is non-unique,
     *  copyForEdit() is called and a new editor is constructed with makeEditor().
     */
    Editor & _edit();

    /**
     *  @Brief Implementation for reserve(); moved here so subclasses can cast the public return value.
     *
     *  _reserve() operates mostly like _edit, but if the capacity is not equal to the current capacity
     *  it will call copyForEdit() with the new capacity even if the editor is unique.
     */
    Editor & _reserve(int capacity);

    /**
     *  @brief Clone the table.
     *
     *  Subclasses should generally implement this with a simple call to the copy constructor,
     *  and shallow-copy as much of their internals as possible.
     */
    virtual Ptr _clone() const = 0;
    
    /**
     *  @brief Called by _edit() when _editor is not unique.
     *
     *  @param[in] capacity    Full size of for new arrays to be allocated.
     *
     *  This should deep-copy all data members in-place.
     *  Subclasses should call the base class implementation, which reallocates
     *  the parameter and weight arrays and throws InvalidParameterException if the
     *  capacity is smaller than the current size.
     */
    virtual void copyForEdit(int capacity);

    /**
     *  @brief Construct a new Editor object that modifies this.
     */
    virtual boost::shared_ptr<Editor> makeEditor() = 0;

private:

    void operator=(SampleTable const &) {} // disabled

    int _size;
    boost::shared_ptr<Editor> _editor;
    lsst::ndarray::Array<double,2,2> _parameters;
    lsst::ndarray::Array<double,1,1> _weights;
};

/**
 *  @brief A separate interface class for mutating operations on SampleTable.
 *
 *  All mutating operations on a subclass SampleTable should be defined as members of
 *  a subclass of Editor.  Editor instances should only be available to users as
 *  references, and should only be constructed by SampleTable::makeEditor.
 *
 *  @todo This should really be an inner class, but SWIG hates those.
 */
class SampleTableEditor : private boost::noncopyable {
public:

    lsst::ndarray::Array<double,1,1> getWeights() {
        return _table->_weights[ndarray::view(0, _table->_size)];
    }

    virtual ~SampleTableEditor() {}
    
protected:

    explicit SampleTableEditor(SampleTable * table) : _table(table) {}

    /**
     *  @brief Append a vector of parameters with associated weights.
     *
     *  Intended to be used by subclasses that implement a public append() method.
     *  Throws LengthErrorException if the capacity of the table is equal to its size.
     */
    void _append(ndarray::Array<double const,1,1> const & parameters, double weight);

private:
    SampleTable * _table;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleTable
