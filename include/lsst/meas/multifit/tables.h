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
#ifndef MEAS_MULTIFIT_tables_h_INCLUDED
#define MEAS_MULTIFIT_tables_h_INCLUDED

#include "lsst/afw/table/Simple.h"
#include "lsst/afw/table/SortedCatalog.h"
#include "lsst/afw/table/BaseColumnView.h"
#include "lsst/afw/table/io/FitsWriter.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/BaseSampler.h"

namespace lsst { namespace meas { namespace multifit {

class ModelFitTable;

/**
 *  @brief Record class used to store galaxy model fitting results
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

    PTR(SampleSet) getSamples() const { return _samples; }
    void setSamples(PTR(SampleSet) samples) { _samples = samples; }

    PTR(afw::detection::Footprint) getFootprint() const { return _footprint; }
    void setFootprint(PTR(afw::detection::Footprint) footprint) { _footprint = footprint; }

protected:

    ModelFitRecord(PTR(ModelFitTable) const & table);

    virtual void _assign(afw::table::BaseRecord const & other);

private:

    PTR(SampleSet) _samples;
    PTR(afw::detection::Footprint) _footprint;
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
     *  @param[in] schema   Schema that defines the fields, offsets, and record size for the table.
     */
    static PTR(ModelFitTable) make(afw::table::Schema const & schema);

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

    ModelFitTable(afw::table::Schema const & schema);

    ModelFitTable(ModelFitTable const & other);

private:

    friend class afw::table::io::FitsWriter;

    template <typename RecordT> friend class afw::table::SortedCatalogT;

     // Return a writer object that knows how to save in FITS format.  See also FitsWriter.
    virtual PTR(afw::table::io::FitsWriter) makeFitsWriter(afw::fits::Fits * fitsfile) const;
};

#ifndef SWIG

typedef afw::table::ColumnViewT<ModelFitRecord> ModelFitColumnView;
typedef afw::table::SortedCatalogT<ModelFitRecord> ModelFitCatalog;
typedef afw::table::SortedCatalogT<ModelFitRecord const> ConstModelFitCatalog;

#endif // !SWIG

}}} // namespace lsst::afw::table

#endif // !MEAS_MULTIFIT_tables_h_INCLUDED
