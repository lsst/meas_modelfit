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

#include <map>
#include <list>

#include "Eigen/LU"

#include "ndarray/eigen.h"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/meas/multifit/parameters.h"

namespace tbl = lsst::afw::table;

namespace lsst { namespace meas { namespace multifit {

namespace {

inline bool isPointNaN(afw::geom::Point2D const & p) {
    return utils::isnan(p.getX()) && utils::isnan(p.getY());
}

class EllipseCoreParameterConverter : public ParameterConverter {
public:

    virtual void apply(
        ndarray::Array<double const,1,1> const & input,
        ndarray::Array<double,1,1> const & output
    ) const {
        _inEllipse->setParameterVector(input.asEigen().segment<3>(0));
        *_outEllipse = *_inEllipse;
        output.asEigen().segment<3>(0) = _outEllipse->getParameterVector();
    }

    virtual double computeJacobian(ndarray::Array<double const,1,1> const & input) const {
        _inEllipse->setParameterVector(input.asEigen().segment<3>(0));
        afw::geom::ellipses::BaseCore::Jacobian j = _outEllipse->dAssign(*_inEllipse);
        return j.determinant();
    }

    EllipseCoreParameterConverter(
        PTR(afw::geom::ellipses::BaseCore) inEllipse,
        PTR(afw::geom::ellipses::BaseCore) outEllipse
    ) : _inEllipse(inEllipse), _outEllipse(outEllipse) {}

private:
    PTR(afw::geom::ellipses::BaseCore) _inEllipse;
    PTR(afw::geom::ellipses::BaseCore) _outEllipse;
};

class EllipseCoreParameterDefinition : public ParameterDefinition {
public:

    virtual PTR(ParameterConverter const) makeConverterTo(ParameterDefinition const & other) const {
        EllipseCoreParameterDefinition const * other2
            = dynamic_cast<EllipseCoreParameterDefinition const *>(&other);
        PTR(ParameterConverter const) result;
        if (!other2) return result;
        result = boost::make_shared<EllipseCoreParameterConverter>(_ellipse, other2->_ellipse);
        return result;
    }

    virtual afw::geom::ellipses::Ellipse makeEllipse(
        ndarray::Array<double const,1,1> const & input
    ) const {
        _ellipse->setParameterVector(input.asEigen().segment<3>(0));
        return afw::geom::ellipses::Ellipse(*_ellipse, _center);
    }

    EllipseCoreParameterDefinition(
        std::string const & name, afw::geom::Point2D const & center
    ) :
        ParameterDefinition(3), _ellipse(afw::geom::ellipses::BaseCore::make(name)), _center(center)
    {}

    virtual bool isPersistable() const { return true; }

protected:

    virtual std::string getPythonModule() const { return "lsst.meas.multifit"; }

    virtual std::string getPersistenceName() const;

    virtual void write(OutputArchiveHandle & handle) const;

    virtual bool _isEqualTo(ParameterDefinition const & other) const {
        if (&other == this) return true;
        EllipseCoreParameterDefinition const & other2
            = static_cast<EllipseCoreParameterDefinition const &>(other);
        return other2._ellipse->getName() == this->_ellipse->getName() &&
            ((isPointNaN(other2._center) || isPointNaN(this->_center)) || other2._center == this->_center);
    }

private:
    PTR(afw::geom::ellipses::BaseCore) _ellipse;
    afw::geom::Point2D _center;
};

class EllipseCoreParameterDefinitionPersistenceKeys : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<std::string> name;
    tbl::Key<tbl::Point<double> > center;

    static EllipseCoreParameterDefinitionPersistenceKeys const & get() {
        static EllipseCoreParameterDefinitionPersistenceKeys const instance;
        return instance;
    }

private:
    EllipseCoreParameterDefinitionPersistenceKeys() :
        schema(),
        name(schema.addField<std::string>("name", "name of the ellipses::BaseCore class", 48)),
        center(schema.addField<tbl::Point<double> >("center",
                                                    "fixed center used when constructing a full ellipse"))
    {
        schema.getCitizen().markPersistent();
    }
};

class EllipseCoreParameterDefinitionFactory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        EllipseCoreParameterDefinitionPersistenceKeys const & keys
            = EllipseCoreParameterDefinitionPersistenceKeys::get();
        LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys.schema);
        tbl::BaseRecord const & record = catalogs.front().front();
        return boost::const_pointer_cast<ParameterDefinition>(
            ParameterDefinition::makeEllipseCoreDefinition(record.get(keys.name), record.get(keys.center))
        );
    }

    explicit EllipseCoreParameterDefinitionFactory(std::string const & name) :
        tbl::io::PersistableFactory(name)
    {}

};


static std::string getEllipseCoreParameterDefinitionPersistenceName() {
    return "EllipseCoreParameterDefinition";
}

// constructor for this instance registers the factor in a singleton in afw::table::io
static EllipseCoreParameterDefinitionFactory registration(getEllipseCoreParameterDefinitionPersistenceName());

std::string EllipseCoreParameterDefinition::getPersistenceName() const {
    return getEllipseCoreParameterDefinitionPersistenceName();
}

void EllipseCoreParameterDefinition::write(OutputArchiveHandle & handle) const {
    EllipseCoreParameterDefinitionPersistenceKeys const & keys
        = EllipseCoreParameterDefinitionPersistenceKeys::get();
    tbl::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    PTR(tbl::BaseRecord) record = catalog.addNew();
    record->set(keys.name, _ellipse->getName());
    record->set(keys.center, _center);
    handle.saveCatalog(catalog);
}

} // anonymous

PTR(ParameterDefinition const) ParameterDefinition::makeEllipseCoreDefinition(
    std::string const & name, afw::geom::Point2D const & center
) {
    return boost::make_shared<EllipseCoreParameterDefinition>(name, center);
}

bool ParameterDefinition::operator==(ParameterDefinition const & other) const {
    if (this->getDim() != other.getDim()) return false;
    if (typeid(*this) != typeid(other)) return false;
    return this->_isEqualTo(other);
}

}}} // namespace lsst::meas::multifit
