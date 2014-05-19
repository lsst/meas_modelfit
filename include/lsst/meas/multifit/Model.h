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

#ifndef LSST_MEAS_MULTIFIT_Model_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Model_h_INCLUDED

#include <vector>

#include "boost/noncopyable.hpp"

#include "lsst/base.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/common.h"

namespace lsst { namespace meas { namespace multifit {

class Model;
class Prior;

typedef std::vector<PTR(Model)> ModelVector;

class Model : private boost::noncopyable {
public:

    enum CenterEnum {
        FIXED_CENTER  = 0x0,
        SINGLE_CENTER = 0x1,
        MULTI_CENTER  = 0x2
    };

    typedef std::vector<std::string> NameVector;
    typedef std::vector<PTR(shapelet::MultiShapeletBasis)> BasisVector;
    typedef std::vector<afw::geom::ellipses::Ellipse> EllipseVector;
    typedef std::vector<afw::geom::ellipses::Ellipse>::iterator EllipseIterator;
    typedef std::vector<afw::geom::ellipses::Ellipse>::const_iterator EllipseConstIterator;

    static PTR(Model) make(BasisVector basisVector, NameVector const & prefixes, CenterEnum center);

    static PTR(Model) make(PTR(shapelet::MultiShapeletBasis) basis, CenterEnum center);

    static PTR(Model) makeGaussian(CenterEnum center, double radius=1.0);

    int getNonlinearDim() const { return _nonlinearNames.size(); }
    int getAmplitudeDim() const { return _amplitudeNames.size(); }
    int getFixedDim() const { return _fixedNames.size(); }
    int getBasisCount() const { return _basisVector.size(); }

    NameVector const & getNonlinearNames() const { return _nonlinearNames; }
    NameVector const & getAmplitudeNames() const { return _amplitudeNames; }
    NameVector const & getFixedNames() const { return _fixedNames; }

    BasisVector const & getBasisVector() const { return _basisVector; }

    shapelet::MultiShapeletFunction makeShapeletFunction(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar const,1,1> const & fixed
    ) const;

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const = 0;

    virtual EllipseVector makeEllipseVector() const = 0;

#ifndef SWIG

    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const = 0;

    virtual void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const = 0;

#endif // !SWIG

    EllipseVector writeEllipses(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & fixed
    ) const;

    void readEllipses(
        EllipseVector const & ellipses,
        ndarray::Array<Scalar,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & fixed
    ) const;

    virtual ~Model() {}

protected:

    Model(
        BasisVector basisVector,
        NameVector nonlinearNames,
        NameVector amplitudeNames,
        NameVector fixedNames
    );

private:
    NameVector _nonlinearNames;
    NameVector _amplitudeNames;
    NameVector _fixedNames;
    BasisVector _basisVector;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Model_h_INCLUDED
