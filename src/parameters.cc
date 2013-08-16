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

#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/parameters.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

typedef std::map<std::string,ParameterDefinition const *> RegistryMap;

RegistryMap & getRegistry() {
    static RegistryMap instance;
    return instance;
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

    static ParameterDefinition const & make(PTR(afw::geom::ellipses::BaseCore) ellipse);

    virtual PTR(ParameterConverter const) makeConverterTo(ParameterDefinition const & other) const {
        EllipseCoreParameterDefinition const * other2
            = dynamic_cast<EllipseCoreParameterDefinition const *>(&other);
        PTR(ParameterConverter const) result;
        if (!other2) return result;
        result = boost::make_shared<EllipseCoreParameterConverter>(_ellipse, other2->_ellipse);
        return result;
    }

    virtual afw::geom::ellipses::Ellipse makeEllipse(
        ndarray::Array<double const,1,1> const & input,
        afw::geom::Point2D const & center=afw::geom::Point2D()
    ) const {
        _ellipse->setParameterVector(input.asEigen().segment<3>(0));
        return afw::geom::ellipses::Ellipse(*_ellipse, center);
    }

    explicit EllipseCoreParameterDefinition(PTR(afw::geom::ellipses::BaseCore) ellipse) :
        ParameterDefinition(ellipse->getName(), 3), _ellipse(ellipse) {}

private:
    PTR(afw::geom::ellipses::BaseCore) _ellipse;
};

ParameterDefinition const & EllipseCoreParameterDefinition::make(
    PTR(afw::geom::ellipses::BaseCore) ellipse
) {
    // Need to keep a static list in order to keep instances alive, since they're constructed on
    // first use, but the main ParameterDefinition registry only stores pointers.
    static std::list<PTR(EllipseCoreParameterDefinition)> instances;
    // TODO: with C++11, use emplace_back() instead of list<PTR(...)>
    instances.push_back(boost::make_shared<EllipseCoreParameterDefinition>(ellipse));
    return *instances.back();
}

} // anonymous

ParameterDefinition const & ParameterDefinition::lookup(std::string const & name_) {
    RegistryMap & registry = getRegistry();
    RegistryMap::const_iterator i = registry.find(name_);
    if (i == registry.end()) {
        try {
            // it should probably throw NotFoundException, but it throws InvalidParameterException instead
            PTR(afw::geom::ellipses::BaseCore) ellipse = afw::geom::ellipses::BaseCore::make(name_);
            return EllipseCoreParameterDefinition::make(ellipse);
        } catch (pex::exceptions::InvalidParameterException &) {}
        throw LSST_EXCEPT(
            pex::exceptions::NotFoundException,
            (boost::format("ParameterDefinition with name '%s' not found in registry")
             % name_).str()
        );
    }
    return *i->second;
}

ParameterDefinition::ParameterDefinition(std::string const & name_, int const size_) :
    name(name_), size(size_)
{
    RegistryMap & registry = getRegistry();
    registry.insert(registry.end(), std::make_pair(name, this));
}

}}} // namespace lsst::meas::multifit
