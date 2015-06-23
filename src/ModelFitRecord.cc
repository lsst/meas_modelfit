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

#include "lsst/afw/table/io/FitsWriter.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"

#include "lsst/meas/modelfit/ModelFitRecord.h"

namespace lsst { namespace meas { namespace modelfit {

//-----------------------------------------------------------------------------------------------------------
//----- Private ModelFitTable/Record classes ---------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// These private derived classes are what you actually get when you do ModelFitTable::make; like the
// private classes in BaseTable.cc, it's more convenient to have an extra set of trivial derived
// classes than to do a lot of friending.

namespace {

class ModelFitTableImpl;

class ModelFitRecordImpl : public ModelFitRecord {
public:

    explicit ModelFitRecordImpl(PTR(ModelFitTable) const & table) : ModelFitRecord(table) {}

};

class ModelFitTableImpl : public ModelFitTable {
public:

    explicit ModelFitTableImpl(
        afw::table::Schema const & schema,
        PTR(afw::table::BaseTable) sampleTable,
        PTR(Interpreter) interpreter
    ) : ModelFitTable(schema, sampleTable, interpreter) {}

    ModelFitTableImpl(ModelFitTableImpl const & other) : ModelFitTable(other) {}

private:

    virtual PTR(afw::table::BaseTable) _clone() const {
        return boost::make_shared<ModelFitTableImpl>(*this);
    }

    virtual PTR(afw::table::BaseRecord) _makeRecord() {
        return boost::make_shared<ModelFitRecordImpl>(getSelf<ModelFitTableImpl>());
    }

};

} // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ModelFitTable/Record member function implementations -----------------------------------------------
//-----------------------------------------------------------------------------------------------------------

ModelFitRecord::ModelFitRecord(PTR(ModelFitTable) const & table) :
    afw::table::SimpleRecord(table), _samples(table->getSampleTable())
{
    if (!table->getSampleTable()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Cannot create ModelFitRecords when the associated SampleTable is not defined "
            "(please call ModelFitTable::setSampleTable)"
        );
    }
}

void ModelFitRecord::_assign(afw::table::BaseRecord const & other) {
    try {
        ModelFitRecord const & s = dynamic_cast<ModelFitRecord const &>(other);
        _footprint = s._footprint;
        _pdf = s._pdf;
        if (s.getSamples().getSchema() == _samples.getSchema()) {
            // We should probably provide the user explicit control over whether to copy the samples
            // (and the pdf and footprint, for that matter), rather than relying on whether it's possible
            // given the schemas.  But that would involve changes to afw::table APIs (probably something
            // like the I/O flags that control how to read/write SourceRecord Footprints), and this works
            // well-enough for now.
            _samples.assign(s.getSamples().begin(), s.getSamples().end(), true);
        }
    } catch (std::bad_cast&) {}
}

PTR(ModelFitTable) ModelFitTable::make(
    afw::table::Schema const & schema,
    PTR(afw::table::BaseTable) sampleTable,
    PTR(Interpreter) interpreter
) {
    if (!checkSchema(schema)) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "Schema for ModelFit must contain at least the keys defined by makeMinimalSchema()."
        );
    }
    return boost::make_shared<ModelFitTableImpl>(schema, sampleTable, interpreter);
}

ModelFitTable::ModelFitTable(
    afw::table::Schema const & schema,
    PTR(afw::table::BaseTable) sampleTable,
    PTR(Interpreter) interpreter
) :
    afw::table::SimpleTable(schema, PTR(afw::table::IdFactory)()),
    _sampleTable(sampleTable),
    _interpreter(interpreter)
{}

ModelFitTable::ModelFitTable(ModelFitTable const & other) :
    afw::table::SimpleTable(other),
    _sampleTable(other._sampleTable->clone()),
    _interpreter(other._interpreter)
{}

PTR(afw::table::io::FitsWriter)
ModelFitTable::makeFitsWriter(afw::fits::Fits * fitsfile, int flags) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicError,
        "FITS persistence is not supported for ModelFitTable"
    );
}

}}} // namespace lsst::meas::modelfit

//-----------------------------------------------------------------------------------------------------------
//----- Explicit instantiation ------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

namespace lsst { namespace afw { namespace table {

template class CatalogT<meas::modelfit::ModelFitRecord>;
template class CatalogT<meas::modelfit::ModelFitRecord const>;

template class SortedCatalogT<meas::modelfit::ModelFitRecord>;
template class SortedCatalogT<meas::modelfit::ModelFitRecord const>;

}}} // namespace lsst::afw::table
