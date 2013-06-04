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

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/BaseTable.h"
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
class ModelFitRecord : public afw::table::BaseRecord {
public:

    typedef ModelFitTable Table;
    typedef afw::table::ColumnViewT<ModelFitRecord> ColumnView;
    typedef afw::table::SortedCatalogT<ModelFitRecord> Catalog;
    typedef afw::table::SortedCatalogT<ModelFitRecord const> ConstCatalog;

    CONST_PTR(ModelFitTable) getTable() const {
        return boost::static_pointer_cast<ModelFitTable const>(afw::table::BaseRecord::getTable());
    }

    afw::table::RecordId getId() const;
    void setId(afw::table::RecordId id);

    afw::geom::ellipses::Quadrupole getShape() const;
    void setShape(afw::geom::ellipses::Quadrupole const & shape);

    afw::geom::Point2D getCentroid() const;
    void setCentroid(afw::geom::Point2D const & centroid);

    double getX() const;
    double getY() const;

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
class ModelFitTable : public afw::table::BaseTable {
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

    /**
     *  @brief Return a minimal schema for ModelFit tables and records.
     *
     *  The returned schema can and generally should be modified further,
     *  but many operations on ModelFitRecords will assume that at least the fields
     *  provided by this routine are present.
     */
    static afw::table::Schema makeMinimalSchema() { return getMinimalSchema().schema; }

    /**
     *  @brief Return true if the given schema is a valid ModelFitTable schema.
     *
     *  This will always be true if the given schema was originally constructed
     *  using makeMinimalSchema(), and will rarely be true otherwise.
     */
    static bool checkSchema(afw::table::Schema const & other) {
        return other.contains(getMinimalSchema().schema);
    }

    /// @brief Key for the unique ID.
    static afw::table::Key<afw::table::RecordId> getIdKey() { return getMinimalSchema().id; }

    /// @brief Key for the centroid
    static afw::table::Key< afw::table::Point<double> > getCentroidKey() {
        return getMinimalSchema().centroid;
    }

    /// @brief Key for the shape
    static afw::table::Key< afw::table::Moments<double> > getShapeKey() {
        return getMinimalSchema().shape;
    }

    /// @copydoc afw::table::BaseTable::clone
    PTR(ModelFitTable) clone() const { return boost::static_pointer_cast<ModelFitTable>(_clone()); }

    /// @copydoc afw::table::BaseTable::makeRecord
    PTR(ModelFitRecord) makeRecord() { return boost::static_pointer_cast<ModelFitRecord>(_makeRecord()); }

    /// @copydoc afw::table::BaseTable::copyRecord
    PTR(ModelFitRecord) copyRecord(afw::table::BaseRecord const & other) {
        return boost::static_pointer_cast<ModelFitRecord>(BaseTable::copyRecord(other));
    }

    /// @copydoc afw::table::BaseTable::copyRecord
    PTR(ModelFitRecord) copyRecord(
        afw::table::BaseRecord const & other,
        afw::table::SchemaMapper const & mapper
    ) {
        return boost::static_pointer_cast<ModelFitRecord>(afw::table::BaseTable::copyRecord(other, mapper));
    }

protected:

    ModelFitTable(afw::table::Schema const & schema);

    ModelFitTable(ModelFitTable const & other);

private:

    // Struct that holds the minimal schema and the special keys we've added to it.
    struct MinimalSchema {
        afw::table::Schema schema;
        afw::table::Key<afw::table::RecordId> id;
        afw::table::Key< afw::table::Point<double> > centroid;
        afw::table::Key< afw::table::Moments<double> > shape;

        MinimalSchema();
    };

    // Return the singleton minimal schema.
    static MinimalSchema & getMinimalSchema();

    friend class afw::table::io::FitsWriter;

    template <typename RecordT> friend class afw::table::SortedCatalogT;

     // Return a writer object that knows how to save in FITS format.  See also FitsWriter.
    virtual PTR(afw::table::io::FitsWriter) makeFitsWriter(afw::fits::Fits * fitsfile, int flags) const;
};

#ifndef SWIG

typedef afw::table::ColumnViewT<ModelFitRecord> ModelFitColumnView;
typedef afw::table::SortedCatalogT<ModelFitRecord> ModelFitCatalog;
typedef afw::table::SortedCatalogT<ModelFitRecord const> ConstModelFitCatalog;

inline afw::table::RecordId ModelFitRecord::getId() const { return get(ModelFitTable::getIdKey()); }
inline void ModelFitRecord::setId(afw::table::RecordId id) { set(ModelFitTable::getIdKey(), id); }

inline afw::geom::ellipses::Quadrupole ModelFitRecord::getShape() const {
    return get(ModelFitTable::getShapeKey());
}
inline void ModelFitRecord::setShape(afw::geom::ellipses::Quadrupole const & shape) {
    return set(ModelFitTable::getShapeKey(), shape);
}

inline afw::geom::Point2D ModelFitRecord::getCentroid() const {
   return get(ModelFitTable::getCentroidKey());
}
inline void ModelFitRecord::setCentroid(afw::geom::Point2D const & centroid) {
    return set(ModelFitTable::getCentroidKey(), centroid);
}

inline double ModelFitRecord::getX() const {
   return get(ModelFitTable::getCentroidKey().getX());
}
inline double ModelFitRecord::getY() const {
   return get(ModelFitTable::getCentroidKey().getY());
}

#endif // !SWIG

}}} // namespace lsst::afw::table

#endif // !MEAS_MULTIFIT_tables_h_INCLUDED
