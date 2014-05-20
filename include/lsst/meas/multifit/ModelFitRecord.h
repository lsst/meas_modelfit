// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
#ifndef MEAS_MULTIFIT_ModelFitRecord_h_INCLUDED
#define MEAS_MULTIFIT_ModelFitRecord_h_INCLUDED

#include "lsst/afw/table/Simple.h"
#include "lsst/afw/table/SortedCatalog.h"
#include "lsst/afw/table/BaseColumnView.h"
#include "lsst/afw/table/io/FitsWriter.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/Mixture.h"
#include "lsst/meas/multifit/Interpreter.h"

namespace lsst { namespace meas { namespace multifit {

class ModelFitTable;

/**
 *  @brief Record class used to store galaxy model fitting results
 *
 *  ModelFitRecord adds a Footprint (indicating the pixels used in the fit,
 *  which is not necessarily the detection Footprint) to each record, and a
 *  Mixture object that provides an analytic approximation to the likelihood
 *  or posterior.  It also joins to a catalog of samples, in a many (samples) to
 *  one (model fit) sense, which can be used to store Monte Carlo or grid
 *  samples from the likelihood or posterior, or the sequence of points taken by
 *  a greedy optimizer.
 */
class ModelFitRecord : public afw::table::SimpleRecord {
public:

    typedef ModelFitTable Table;
    typedef afw::table::ColumnViewT<ModelFitRecord> ColumnView;
    typedef afw::table::SortedCatalogT<ModelFitRecord> Catalog;
    typedef afw::table::SortedCatalogT<ModelFitRecord const> ConstCatalog;

    CONST_PTR(ModelFitTable) getTable() const {
        return boost::static_pointer_cast<ModelFitTable const>(afw::table::BaseRecord::getTable());
    }

    PTR(Interpreter) getInterpreter() const;

    afw::table::BaseCatalog const & getSamples() const { return _samples; }
    afw::table::BaseCatalog & getSamples() { return _samples; }

    PTR(afw::detection::Footprint) getFootprint() const { return _footprint; }
    void setFootprint(PTR(afw::detection::Footprint) footprint) { _footprint = footprint; }

    PTR(Mixture) getPdf() const { return _pdf; }
    void setPdf(PTR(Mixture) pdf) { _pdf = pdf; }

protected:

    ModelFitRecord(PTR(ModelFitTable) const & table);

    virtual void _assign(afw::table::BaseRecord const & other);

private:

    afw::table::BaseCatalog _samples;
    PTR(afw::detection::Footprint) _footprint;
    PTR(Mixture) _pdf;
};

/**
 *  @brief Table class used to store galaxy model fitting results
 *
 *  @copydetails ModelFitRecord
 */
class ModelFitTable : public afw::table::SimpleTable {
public:

    typedef ModelFitRecord Record;
    typedef afw::table::ColumnViewT<ModelFitRecord> ColumnView;
    typedef afw::table::SortedCatalogT<Record> Catalog;
    typedef afw::table::SortedCatalogT<Record const> ConstCatalog;

    /**
     *  @brief Construct a new table.
     *
     *  @param[in] schema   Schema that defines the fields, offsets, and record size for the main table.
     *  @param[in] sampleTable   A table object for the samples associated with the fitting results
     *                           in a many (samples) to one (model fits) relationship.
     *  @param[in] interpreter   An instance of Interpreter that can be used to interpret the
     *                           Sample catalog and Pdf Mixture.
     */
    static PTR(ModelFitTable) make(
        afw::table::Schema const & schema,
        PTR(afw::table::BaseTable) sampleTable = PTR(afw::table::BaseTable)(),
        PTR(Interpreter) interpreter = PTR(Interpreter)()
    );

    /// Return the table object used to allocate records in the related sample catalogs.
    PTR(afw::table::BaseTable) getSampleTable() const { return _sampleTable; }

    /// Set the table object used to allocate records in the related sample catalogs.
    void setSampleTable(PTR(afw::table::BaseTable) sampleTable) { _sampleTable = sampleTable; }

    /// Return an object that can be used to interpret the attached Samples and Pdf
    PTR(Interpreter) getInterpreter() const { return _interpreter; }

    /// Set the object that can be used to interpret the attached Samples and Pdf
    void setInterpreter(PTR(Interpreter) interpreter) { _interpreter = interpreter; }

    /// @copydoc afw::table::BaseTable::clone
    PTR(ModelFitTable) clone() const { return boost::static_pointer_cast<ModelFitTable>(_clone()); }

    /// @copydoc afw::table::BaseTable::makeRecord
    PTR(ModelFitRecord) makeRecord() { return boost::static_pointer_cast<ModelFitRecord>(_makeRecord()); }

    /// @copydoc afw::table::BaseTable::copyRecord
    PTR(ModelFitRecord) copyRecord(afw::table::BaseRecord const & other) {
        return boost::static_pointer_cast<ModelFitRecord>(SimpleTable::copyRecord(other));
    }

    /// @copydoc afw::table::BaseTable::copyRecord
    PTR(ModelFitRecord) copyRecord(
        afw::table::BaseRecord const & other,
        afw::table::SchemaMapper const & mapper
    ) {
        return boost::static_pointer_cast<ModelFitRecord>(afw::table::SimpleTable::copyRecord(other, mapper));
    }

protected:

    ModelFitTable(
        afw::table::Schema const & schema,
        PTR(afw::table::BaseTable) sampleTable,
        PTR(Interpreter) interpreter
    );

    ModelFitTable(ModelFitTable const & other);

private:

    friend class afw::table::io::FitsWriter;

    template <typename RecordT> friend class afw::table::SortedCatalogT;

     // Return a writer object that knows how to save in FITS format.  See also FitsWriter.
    virtual PTR(afw::table::io::FitsWriter) makeFitsWriter(afw::fits::Fits * fitsfile, int flags) const;

    PTR(afw::table::BaseTable) _sampleTable;
    PTR(Interpreter) _interpreter;
};

#ifndef SWIG

inline PTR(Interpreter) ModelFitRecord::getInterpreter() const { return getTable()->getInterpreter(); }

typedef afw::table::ColumnViewT<ModelFitRecord> ModelFitColumnView;
typedef afw::table::SortedCatalogT<ModelFitRecord> ModelFitCatalog;
typedef afw::table::SortedCatalogT<ModelFitRecord const> ConstModelFitCatalog;

#endif // !SWIG

}}} // namespace lsst::afw::table

#endif // !MEAS_MULTIFIT_ModelFitRecord_h_INCLUDED
