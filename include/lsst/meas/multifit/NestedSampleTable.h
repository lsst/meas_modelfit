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

#ifndef LSST_MEAS_MULTIFIT_NestedSampleTable
#define LSST_MEAS_MULTIFIT_NestedSampleTable

#include "lsst/meas/multifit/SampleTable.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  An intermediate base class SampleTable that adds a nested Gaussian at each sample point.
 *
 *  Each nested Gaussian is represented by an amplitude, a mean vector, and a matrix that
 *  may be the covariance matrix, its inverse, or the lower Cholesky factor of either.
 */
class NestedSampleTable : public SampleTable {
public:

    typedef boost::shared_ptr<NestedSampleTable> Ptr;
    typedef boost::shared_ptr<NestedSampleTable const> ConstPtr;

    enum NestedMatrixType {
        COVARIANCE,     ///< Standard covariance matrix.
        FISHER,         ///< Fisher information matrix == inverse of covariance.
        COVARIANCE_LLT, ///< Lower-triangular Cholesky factor of the covariance.
        FISHER_LLT      ///< Lower-triangular Cholesky factor of the Fisher matrix.
    };

    /**
     *  @brief The vector of nested amplitudes.
     */
    lsst::ndarray::Array<double const,1,1> getNestedAmplitudes() const {
        return _nestedAmplitudes[ndarray::view(0, getSize())];
    }

    /**
     *  @brief The (samples)x(nested dim.) matrix of nested mean vectors.
     */
    lsst::ndarray::Array<double const,2,2> getNestedMeans() const {
        return _nestedMeans[ndarray::view(0, getSize())];
    }

    /**
     *  @brief The (samples)x(nested dim.)x(nested dim.) tensor of nested matrices.
     */
    lsst::ndarray::Array<double const,3,3> getNestedMatrices() const {
        return _nestedMatrices[ndarray::view(0, getSize())];
    }

    /// @brief Return how the nested matrices should be interpreted.
    NestedMatrixType getNestedMatrixType() const { return _nestedMatrixType; }

    /// @brief Return the number of parameters for the nested Gaussian distributions.
    int getNestedDimensionality() const { return _nestedMeans.getSize<1>(); }

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

        lsst::ndarray::Array<double,1,1> const & getNestedAmplitudes() {
            return getTable()._nestedAmplitudes;
        }

        lsst::ndarray::Array<double,2,2> const & getNestedMeans() {
            return getTable()._nestedMeans;
        }

        lsst::ndarray::Array<double,3,3> const & getNestedMatrices() {
            return getTable()._nestedMatrices;
        }

        NestedSampleTable & getTable() {
            return static_cast<NestedSampleTable&>(SampleTable::Editor::getTable());
        }
    };
#endif

    /// @brief Construct with given capacity and dimensionalities.
    NestedSampleTable(int capacity, int dimensionality, int nestedDimensionality,
                      NestedMatrixType nestedMatrixType);

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
    NestedMatrixType _nestedMatrixType;
    lsst::ndarray::Array<double,1,1> _nestedAmplitudes;
    lsst::ndarray::Array<double,2,2> _nestedMeans;
    lsst::ndarray::Array<double,3,3> _nestedMatrices;
};


}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleTable
