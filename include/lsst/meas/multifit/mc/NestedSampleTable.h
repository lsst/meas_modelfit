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

#ifndef LSST_MEAS_MULTIFIT_MC_NestedSampleTable
#define LSST_MEAS_MULTIFIT_MC_NestedSampleTable

#include "lsst/meas/multifit/mc/SampleTable.h"

namespace lsst { namespace meas { namespace multifit { namespace mc {

/**
 *  @brief An intermediate base class SampleTable that adds a nested table of weighted coefficient vectors
 *         at each parameter point.
 */
class NestedSampleTable : public SampleTable {
public:

    typedef boost::shared_ptr<NestedSampleTable> Ptr;
    typedef boost::shared_ptr<NestedSampleTable const> ConstPtr;

    /**
     *  @brief The (size)x(nested size) array of nested weights for each coefficient vector.
     */
    lsst::ndarray::Array<double const,2,2> getNestedWeights() const {
        return _nestedWeights[ndarray::view(0, getSize())];
    }

    /**
     *  @brief The (size)x(nested size)(coefficient count) array of nested coefficient values.
     */
    lsst::ndarray::Array<Pixel const,3,3> getCoefficients() const {
        return _coefficients[ndarray::view(0, getSize())];
    }

    /// @brief Return the number of nested coefficient samples for each parameter sample.
    int getNestedSize() const { return _nestedWeights.getSize<1>(); }

    /// @brief Return the number of coefficients.
    int getCoefficientCount() const { return _coefficients.getSize<2>(); }

    /**
     *  @brief Copy the table (shallow with copy-on-write).
     *
     *  Subclasses should reimplement, calling _clone() and casting the returned pointer
     *  to the appropriate subclass of SampleTable.
     */
    Ptr clone() const { return boost::static_pointer_cast<NestedSampleTable>(_clone()); }

protected:

#ifndef SWIG
    class Editor : public SampleTable::Editor {
    protected:

        explicit Editor(NestedSampleTable * table) : SampleTable::Editor(table) {}

        ndarray::Array<double,2,2> const & getNestedWeights() {
            return getTable()._nestedWeights;
        }

        ndarray::Array<Pixel,3,3> const & getCoefficients() {
            return getTable()._coefficients;
        }

        NestedSampleTable & getTable() {
            return static_cast<NestedSampleTable&>(SampleTable::Editor::getTable());
        }
    };
#endif

    /// @brief Construct with given capacity, nested sample size, parameter count, and coefficient count.
    NestedSampleTable(int capacity, int nestedSize, int parameterCount, int coefficientCount);

    /// @brief Copy constructor.
    NestedSampleTable(NestedSampleTable const & other);

    /// @brief Subset copy constructor.
    NestedSampleTable(NestedSampleTable const & other, int start, int stop);
    
    /**
     *  @brief Called by _edit() when _editor is not unique.
     *
     *  @param[in] capacity    Full size of for new arrays to be allocated.
     *
     *  This should deep-copy all data members in-place.
     *  Subclasses should call their immediate base class implementation.
     */
    virtual void copyForEdit(int capacity);

private:
    ndarray::Array<double,2,2> _nestedWeights;
    ndarray::Array<Pixel,3,3> _coefficients;
};


}}}} // namespace lsst::meas::multifit::mc

#endif // !LSST_MEAS_MULTIFIT_MC_SampleTable
