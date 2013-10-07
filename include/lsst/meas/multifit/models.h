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

#ifndef LSST_MEAS_MULTIFIT_models_h_INCLUDED
#define LSST_MEAS_MULTIFIT_models_h_INCLUDED

#include <vector>

#include "boost/noncopyable.hpp"

#include "lsst/base.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/shapelet/MultiShapeletBasis.h"

namespace lsst { namespace meas { namespace multifit {

class Model;
class Prior;

typedef std::vector<PTR(Model)> ModelVector;

class Model : private boost::noncopyable {
public:

    typedef std::vector<PTR(shapelet::MultiShapeletBasis)> BasisVector;
    typedef std::vector<afw::geom::ellipses::Ellipse> EllipseVector;
    typedef std::vector<afw::geom::ellipses::Ellipse>::iterator EllipseIterator;
    typedef std::vector<afw::geom::ellipses::Ellipse>::const_iterator EllipseConstIterator;

    static PTR(Model) makeFixedCenter(BasisVector basisVector, afw::geom::Point2D const & center);
    static PTR(Model) makeSingleCenter(BasisVector basisVector);
    static PTR(Model) makeMultiCenter(BasisVector basisVector);

    int getParameterDim() const { return _parameterDim; }
    int getAmplitudeDim() const { return _amplitudeDim; }
    int getBasisCount() const { return _basisVector.size(); }

    BasisVector const & getBasisVector() const { return _basisVector; }

    shapelet::MultiShapeletFunction makeShapeletFunction(
        ndarray::Array<double const,1,1> const & parameters,
        ndarray::Array<double const,1,1> const & amplitudes
    ) const;

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const = 0;

    virtual EllipseVector makeEllipseVector() const = 0;

    virtual void writeEllipses(double const * parameterIter, EllipseIterator ellipseIter) const = 0;

    virtual void readEllipses(EllipseConstIterator ellipseIter, double * parameterIter) const = 0;

    virtual ~Model() {}

protected:

    Model(BasisVector basisVector, int parameterDim);

private:
    int _parameterDim;
    int _amplitudeDim;
    BasisVector _basisVector;
};

class MultiModel : public Model {
public:

    explicit MultiModel(ModelVector components);

    ModelVector const & getComponents() const { return _components; }

    virtual EllipseVector makeEllipseVector() const;

    virtual void writeEllipses(double const * parameterIter, EllipseIterator ellipseIter) const;

    virtual void readEllipses(EllipseConstIterator ellipseIter, double * parameterIter) const;

private:
    ModelVector _components;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_models_h_INCLUDED
