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
#include "lsst/meas/multifit/constants.h"
#include "boost/shared_ptr.hpp"
#include "boost/noncopyable.hpp"

namespace lsst { namespace meas { namespace multifit {

/**
 *  A base class for tables of weighted samples of points in parameter space.
 */
template <typename T>
class SampleTable {
public:

    typedef boost::shared_ptr<SampleTable> Ptr;
    typedef boost::shared_ptr<SampleTable const> ConstPtr;

    /**
     *  @brief The (size)x(parameter count) array of parameter values.
     */
    lsst::ndarray::Array<T const,2,2> getParameters() const {
        return _parameters[ndarray::view(0, _size)];
    }

    /**
     *  @brief The weights corresponding to each sample.
     */
    lsst::ndarray::Array<T const,1,1> getWeights() const {
        return _weights[ndarray::view(0, _size)];
    }

    /// @brief Return the number of parameters.
    int getParameterCount() const { return _parameters.getSize<1>(); }

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

    virtual ~SampleTable() {}

protected:

#ifndef SWIG
    /**
     *  @brief A separate interface class for mutating operations on SampleTable.
     *
     *  All mutating operations on a subclass SampleTable should be defined as members of
     *  a subclass of Editor.  Editor instances should only be available to users as
     *  references, and should only be constructed by SampleTable::makeEditor.
     */
    class Editor : private boost::noncopyable {
    public:

        typedef boost::shared_ptr<Editor> Ptr;

        virtual ~Editor() {}
    
    protected:

        explicit Editor(SampleTable * table) : _table(table) {}

        ndarray::Array<T,1,1> const & getWeights() {
            return getTable()._weights;
        }

        ndarray::Array<T,2,2> const & getParameters() {
            return getTable()._parameters;
        }

        int & getSize() { return getTable()._size; }

        SampleTable & getTable() { return *_table; }

    private:
        SampleTable * _table;
    };
#endif

    /// @brief Construct with given capacity and parameter count.
    SampleTable(int capacity, int parameterCount);

    /// @brief Copy constructor.
    SampleTable(SampleTable const & other);

    /// @brief Subset copy constructor.
    SampleTable(SampleTable const & other, int start, int stop);

    /**
     *  @brief Invoke copy-on-write.
     *
     *  The copy-on-write mechanism works by checking the _editor data member.  If it is empty (as
     *  it will be for a newly-constructed, non-copied table) makeEditor() is called to construct
     *  a new Editor.  If it is unique, the existing Editor is returned.  If it is non-unique
     *  or the given capacity is not equal to the current capacity, copyForEdit() is called and
     *  a new editor is constructed with makeEditor().
     */
    Editor & _edit(int capacity);

    Editor & _edit() { return _edit(getCapacity()); }

    /**
     *  @brief Clone the table.
     *
     *  Subclasses should generally implement this with a simple call to the copy constructor,
     *  and shallow-copy as much of their internals as possible.
     */
    virtual Ptr _clone() const = 0;
    
    /**
     *  @brief Called by _edit() and _reserve() when _editor is not unique.
     *
     *  @param[in] capacity    Full size for new arrays to be allocated.
     *
     *  This should deep-copy all data members in-place.
     *  Subclasses should call their immediate base class implementation.
     *  This implementation reallocates the parameter and weight arrays and throws
     *  InvalidParameterException if the capacity is smaller than the current size.
     */
    virtual void copyForEdit(int capacity);

    /**
     *  @brief Construct a new Editor object that modifies this.
     */
    virtual Editor::Ptr makeEditor() = 0;

    /// @brief Helper function to aid in implementing copyForEdit.
    template <typename T, int N, int C>
    void copyArrayForEdit(ndarray::Array<T,N,C> & array, int capacity) const {
        ndarray::Array<T,N,C> newArray(
            ndarray::allocate(ndarray::concatenate(capacity, array.getShape().template first<N-1>()))
        );
        newArray[ndarray::view(0, _size)] = array[ndarray::view(0, _size)];
        array = newArray;
    }

private:

    void operator=(SampleTable const &) {} // disabled

    int _size;
    Editor::Ptr _editor;
    ndarray::Array<T,2,2> _parameters;
    ndarray::Array<T,1,1> _weights;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleTable
