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

#include "Eigen/LU"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/integrals.h"

namespace tbl = lsst::afw::table;

namespace lsst { namespace meas { namespace multifit {

samples::Scalar FlatPrior::apply(LogGaussian const & likelihood, samples::Vector const & parameters) const {
    return integrateGaussian(likelihood.grad, likelihood.fisher)
        + std::log(_maxRadius * _maxEllipticity * _maxEllipticity * 2 * M_PI);
}

FlatPrior::FlatPrior(double maxRadius, double maxEllipticity) :
    Prior(ParameterDefinition::lookup("SeparableReducedShearTraceRadius")),
    _maxRadius(maxRadius),
    _maxEllipticity(maxEllipticity)
{}

namespace {

class FlatPriorPersistenceKeys : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<double> maxRadius;
    tbl::Key<double> maxEllipticity;

    static FlatPriorPersistenceKeys const & get() {
        static FlatPriorPersistenceKeys const instance;
        return instance;
    }
private:
    FlatPriorPersistenceKeys() :
        schema(),
        maxRadius(schema.addField<double>("radius.max", "maximum allowed radius")),
        maxEllipticity(schema.addField<double>("ellipticty.max", "maximum allowed ellipticity"))
    {
        schema.getCitizen().markPersistent();
    }
};

class FlatPriorFactory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        FlatPriorPersistenceKeys const & keys = FlatPriorPersistenceKeys::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys.schema);
        tbl::BaseRecord const & record = catalogs.front().front();
        return boost::make_shared<FlatPrior>(record.get(keys.maxRadius), record.get(keys.maxEllipticity));
    }

    explicit FlatPriorFactory(std::string const & name) : tbl::io::PersistableFactory(name) {}

};

std::string getFlatPriorPersistenceName() { return "FlatPrior"; }

FlatPriorFactory registration(getFlatPriorPersistenceName());

} // anonymous

std::string FlatPrior::getPersistenceName() const { return getFlatPriorPersistenceName(); }

void FlatPrior::write(OutputArchiveHandle & handle) const {
    FlatPriorPersistenceKeys const & keys = FlatPriorPersistenceKeys::get();
    tbl::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    PTR(tbl::BaseRecord) record = catalog.addNew();
    record->set(keys.maxRadius, _maxRadius);
    record->set(keys.maxEllipticity, _maxEllipticity);
    handle.saveCatalog(catalog);
}

}}} // namespace lsst::meas::multifit
