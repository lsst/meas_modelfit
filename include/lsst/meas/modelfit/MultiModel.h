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

#ifndef LSST_MEAS_MODELFIT_MultiModel_h_INCLUDED
#define LSST_MEAS_MODELFIT_MultiModel_h_INCLUDED

#include "lsst/meas/modelfit/Model.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  @brief A concrete Model class that simply concatenates several other Models
 */
class MultiModel : public Model {
public:

    /**
     *  @brief Construct a new MultiModel
     *
     *  @param[in] components       A vector of other Models to combine
     *  @param[in] prefixes         A vector of name prefixes used to construct parameter names,
     *                              one for each element in components.
     */
    explicit MultiModel(ModelVector components, NameVector const & prefixes);

    /// Return the vector of constituent models
    ModelVector const & getComponents() const { return _components; }

    /// @copydoc Model::adaptPrior
    std::shared_ptr<Prior> adaptPrior(std::shared_ptr<Prior> prior) const override;

    /// @copydoc Model::makeEllipseVector
    EllipseVector makeEllipseVector() const override;

    /// @copydoc Model::writeEllipses
    void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const override;

    /// @copydoc Model::readEllipses
    void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const override;

private:
    ModelVector _components;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_MultiModel_h_INCLUDED
