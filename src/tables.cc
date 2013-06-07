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

#include "lsst/meas/multifit/tables.h"

namespace lsst { namespace meas { namespace multifit {

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

    explicit ModelFitTableImpl(afw::table::Schema const & schema) : ModelFitTable(schema) {}

    ModelFitTableImpl(ModelFitTableImpl const & other) : ModelFitTable(other) {}

private:

    virtual PTR(afw::table::BaseTable) _clone() const {
        return boost::make_shared<ModelFitTableImpl>(*this);
    }

    virtual PTR(afw::table::BaseRecord) _makeRecord() {
        return boost::make_shared<ModelFitRecordImpl>(getSelf<ModelFitTableImpl>());
    }

};

// Schema prepended when saving a ModelFit table
struct PersistenceSchema : private boost::noncopyable {
    afw::table::Schema schema;
    afw::table::Key<int> samples;
    afw::table::Key<int> footprint;

    static PersistenceSchema const & get() {
        static PersistenceSchema const instance;
        return instance;
    }

    // Create a SchemaMapper that maps a ModelFitRecord to a BaseRecord with ID for SampleSet
    afw::table::SchemaMapper makeWriteMapper(afw::table::Schema const & inputSchema) const {
        std::vector<afw::table::Schema> inSchemas;
        inSchemas.push_back(PersistenceSchema::get().schema);
        inSchemas.push_back(inputSchema);
        return afw::table::SchemaMapper::join(inSchemas).back(); // don't need front; it's an identity mapper
    }

    // Create a SchemaMapper that maps a BaseRecord with ID for SampleSet to a ModelFitRecord
    afw::table::SchemaMapper makeReadMapper(afw::table::Schema const & inputSchema) const {
        return afw::table::SchemaMapper::removeMinimalSchema(inputSchema, schema);
    }

    // Convert a ModelFitRecord to a BaseRecord with IDs for Psf and Wcs.
    template <typename OutputArchiveIsh>
    void writeRecord(
        ModelFitRecord const & input, afw::table::BaseRecord & output,
        afw::table::SchemaMapper const & mapper, OutputArchiveIsh & archive
    ) const {
        output.assign(input, mapper);
        output.set(samples, archive.put(input.getSamples()));
        output.set(footprint, archive.put(input.getFootprint()));
    }

    void readRecord(
        afw::table::BaseRecord const & input, ModelFitRecord & output,
        afw::table::SchemaMapper const & mapper, afw::table::io::InputArchive const & archive
    ) const {
        output.assign(input, mapper);
        output.setSamples(archive.get<SampleSet>(input.get(samples)));
        output.setFootprint(archive.get<afw::detection::Footprint>(input.get(footprint)));
    }

private:
    PersistenceSchema() :
        schema(),
        samples(schema.addField<int>("samples", "archive ID for SampleSet object")),
        footprint(schema.addField<int>("samplesfootpring", "archive ID for Footprint object"))
    {
        schema.getCitizen().markPersistent();
    }
};

} // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ModelFitFitsWriter ---------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// A custom FitsWriter for ModelFit - this just sets the AFW_TYPE key to MODELFIT, which should ensure
// we use ModelFitFitsReader to read it.

namespace {

class ModelFitFitsWriter : public afw::table::io::FitsWriter {
public:

    ModelFitFitsWriter(Fits * fits, int flags) : afw::table::io::FitsWriter(fits, flags), _archive() {}

protected:

    virtual void _writeTable(CONST_PTR(afw::table::BaseTable) const & table, std::size_t nRows);

    virtual void _writeRecord(afw::table::BaseRecord const & r);

    virtual void _finish() { _archive.writeFits(*_fits); }

    afw::table::io::OutputArchive _archive;
    PTR(afw::table::BaseRecord) _record;
    afw::table::SchemaMapper _mapper;
};

void ModelFitFitsWriter::_writeTable(CONST_PTR(afw::table::BaseTable) const & t, std::size_t nRows) {
    CONST_PTR(ModelFitTable) inTable = boost::dynamic_pointer_cast<ModelFitTable const>(t);
    if (!inTable) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot use a ModelFitFitsWriter on a non-ModelFit table."
        );
    }
    _mapper = PersistenceSchema::get().makeWriteMapper(inTable->getSchema());
    PTR(afw::table::BaseTable) outTable = afw::table::BaseTable::make(_mapper.getOutputSchema());
    afw::table::io::FitsWriter::_writeTable(outTable, nRows);
    _fits->writeKey("AFW_TYPE", "MODELFIT", "Tells lsst::afw to load this as a ModelFit table.");
    _record = outTable->makeRecord();
}

void ModelFitFitsWriter::_writeRecord(afw::table::BaseRecord const & r) {
    ModelFitRecord const & record = static_cast<ModelFitRecord const &>(r);
    PersistenceSchema::get().writeRecord(record, *_record, _mapper, _archive);
    afw::table::io::FitsWriter::_writeRecord(*_record);
}

} // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ModelFitFitsReader ---------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// A custom FitsReader for ModelFitTable/Record - this gets registered with name MODELFIT, so it should get
// used whenever we read a table with AFW_TYPE set to that value.

namespace {

class ModelFitFitsReader : public afw::table::io::FitsReader {
public:

    explicit ModelFitFitsReader(Fits * fits, PTR(afw::table::io::InputArchive) archive, int flags) :
        afw::table::io::FitsReader(fits, archive, flags), _archive(archive)
    {
        if (!_archive) {
            int oldHdu = _fits->getHdu();
            _fits->setHdu(oldHdu + 1);
            _archive.reset(new afw::table::io::InputArchive(afw::table::io::InputArchive::readFits(*_fits)));
            _fits->setHdu(oldHdu);
        }
    }

protected:

    virtual PTR(afw::table::BaseTable) _readTable();

    virtual PTR(afw::table::BaseRecord) _readRecord(PTR(afw::table::BaseTable) const & table);

    PTR(afw::table::BaseTable) _inTable;
    PTR(afw::table::io::InputArchive) _archive;
    afw::table::SchemaMapper _mapper;
};

PTR(afw::table::BaseTable) ModelFitFitsReader::_readTable() {
    PTR(daf::base::PropertyList) metadata = boost::make_shared<daf::base::PropertyList>();
    _fits->readMetadata(*metadata, true);
    afw::table::Schema schema(*metadata, true);
    _inTable = afw::table::BaseTable::make(schema);
    _mapper = PersistenceSchema::get().makeReadMapper(schema);
    PTR(ModelFitTable) table = ModelFitTable::make(_mapper.getOutputSchema());
    _startRecords(*table);
    if (metadata->exists("AFW_TYPE")) metadata->remove("AFW_TYPE");
    table->setMetadata(metadata);
    return table;
}

PTR(afw::table::BaseRecord) ModelFitFitsReader::_readRecord(PTR(afw::table::BaseTable) const & t) {
    PTR(ModelFitRecord) record;
    PTR(ModelFitTable) table = boost::static_pointer_cast<ModelFitTable>(t);
    PTR(afw::table::BaseRecord) inRecord = afw::table::io::FitsReader::_readRecord(_inTable);
    if (inRecord) {
        record = table->makeRecord();
        PersistenceSchema::get().readRecord(*inRecord, *record, _mapper, *_archive);
    }
    return record;
}

// registers the reader so FitsReader::make can use it.
static afw::table::io::FitsReader::FactoryT<ModelFitFitsReader> referenceFitsReaderFactory("MODELFIT");

} // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ModelFitTable/Record member function implementations -----------------------------------------------
//-----------------------------------------------------------------------------------------------------------

ModelFitRecord::ModelFitRecord(PTR(ModelFitTable) const & table) : afw::table::SimpleRecord(table) {}

void ModelFitRecord::_assign(afw::table::BaseRecord const & other) {
    try {
        ModelFitRecord const & s = dynamic_cast<ModelFitRecord const &>(other);
        _samples = s._samples;
        _footprint = s._footprint;
    } catch (std::bad_cast&) {}
}

PTR(ModelFitTable) ModelFitTable::make(afw::table::Schema const & schema) {
    if (!checkSchema(schema)) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "Schema for ModelFit must contain at least the keys defined by makeMinimalSchema()."
        );
    }
    return boost::make_shared<ModelFitTableImpl>(schema);
}

ModelFitTable::ModelFitTable(afw::table::Schema const & schema) :
    afw::table::SimpleTable(schema, PTR(afw::table::IdFactory)()) {}

ModelFitTable::ModelFitTable(ModelFitTable const & other) :
    afw::table::SimpleTable(other) {}

PTR(afw::table::io::FitsWriter)
ModelFitTable::makeFitsWriter(afw::fits::Fits * fitsfile, int flags) const {
    return boost::make_shared<ModelFitFitsWriter>(fitsfile, flags);
}

}}} // namespace lsst::meas::multifit

//-----------------------------------------------------------------------------------------------------------
//----- Explicit instantiation ------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

namespace lsst { namespace afw { namespace table {

template class CatalogT<meas::multifit::ModelFitRecord>;
template class CatalogT<meas::multifit::ModelFitRecord const>;

template class SortedCatalogT<meas::multifit::ModelFitRecord>;
template class SortedCatalogT<meas::multifit::ModelFitRecord const>;

}}} // namespace lsst::afw::table
