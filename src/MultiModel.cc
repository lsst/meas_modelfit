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

#include "lsst/pex/exceptions.h"
#include "lsst/meas/modelfit/MultiModel.h"
#include "lsst/meas/modelfit/Prior.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

static Model::BasisVector concatenateBasisVectors(ModelVector const & components) {
    Model::BasisVector r;
    for (ModelVector::const_iterator i = components.begin(); i != components.end(); ++i) {
        r.insert(r.end(), (**i).getBasisVector().begin(), (**i).getBasisVector().end());
    }
    return r;
}

typedef Model::NameVector const & (Model::*ModelNameGetter)() const;

static Model::NameVector concatenateNameVectors(
    ModelVector const & components, Model::NameVector const & prefixes, ModelNameGetter getter
) {
    LSST_THROW_IF_NE(
        components.size(), prefixes.size(),
        pex::exceptions::LengthError,
        "Number of model components (%d) does not match number of prefixes (%d)"
    );
    Model::NameVector r;
    for (std::size_t i = 0, ni = components.size(); i < ni; ++i) {
        Model::NameVector const & componentNames = ((*components[i]).*getter)();
        for (std::size_t j = 0, nj = componentNames.size(); j < nj; ++j) {
            r.push_back(prefixes[i] + componentNames[j]);
        }
    }
    return r;
}

} // anonymous

MultiModel::MultiModel(ModelVector components, NameVector const & prefixes) :
    Model(
        concatenateBasisVectors(components),
        concatenateNameVectors(components, prefixes, &Model::getNonlinearNames),
        concatenateNameVectors(components, prefixes, &Model::getAmplitudeNames),
        concatenateNameVectors(components, prefixes, &Model::getFixedNames)
    ),
    _components(components)
{}

std::shared_ptr<Prior> MultiModel::adaptPrior(std::shared_ptr<Prior> prior) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicError,
        "adaptPrior not implemented for MultiModel"
    );
}

Model::EllipseVector MultiModel::makeEllipseVector() const {
    EllipseVector r;
    for (ModelVector::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        EllipseVector c = (**i).makeEllipseVector();
        r.insert(r.end(), c.begin(), c.end());
    }
    return r;
}

void MultiModel::writeEllipses(
    Scalar const * nonlinearIter, Scalar const * fixedIter,
    EllipseIterator ellipseIter
) const {
    for (ModelVector::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        (**i).writeEllipses(nonlinearIter, fixedIter, ellipseIter);
        nonlinearIter += (**i).getNonlinearDim();
        fixedIter += (**i).getFixedDim();
        ellipseIter += (**i).getBasisCount();
    }
}

void MultiModel::readEllipses(
    EllipseConstIterator ellipseIter,
    Scalar * nonlinearIter, Scalar * fixedIter
) const {
    for (ModelVector::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        (**i).readEllipses(ellipseIter, nonlinearIter, fixedIter);
        nonlinearIter += (**i).getNonlinearDim();
        fixedIter += (**i).getFixedDim();
        ellipseIter += (**i).getBasisCount();
    }
}

}}} // namespace lsst::meas::modelfit
