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

#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

static afw::geom::ellipses::SeparableConformalShearLogTraceRadius extendedSourceEllipseCore;
static afw::geom::ellipses::Quadrupole pointSourceEllipseCore(0.0, 0.0, 0.0);

int countParameters(Model::BasisVector const & basisVector, int ellipseFactor, int pointFactor) {
    int r = 0;
    for (Model::BasisVector::const_iterator i = basisVector.begin(); i != basisVector.end(); ++i) {
        if (*i) {
            r += ellipseFactor + pointFactor;
        } else {
            r += pointFactor;
        }
    }
    return r;
}

Model::EllipseVector makeEllipseVectorImpl(
    Model::BasisVector const & basisVector,
    afw::geom::Point2D const & center
) {
    Model::EllipseVector r;
    r.reserve(basisVector.size());
    for (int i = 0, n = basisVector.size(); i < n; ++i) {
        if (basisVector[i]) {
            r.push_back(afw::geom::ellipses::Ellipse(extendedSourceEllipseCore, center));
        } else {
            r.push_back(afw::geom::ellipses::Ellipse(pointSourceEllipseCore, center));
        }
    }
    return r;
}

} // anonymous

// ========== Model =========================================================================================

shapelet::MultiShapeletFunction Model::makeShapeletFunction(
    ndarray::Array<double const,1,1> const & parameters,
    ndarray::Array<double const,1,1> const & amplitudes
) const {
    EllipseVector ellipses = makeEllipseVector();
    writeEllipses(parameters.begin(), ellipses.begin());
    shapelet::MultiShapeletFunction r;
    int c = 0;
    for (int i = 0, n = getBasisCount(); i < n; ++i) {
        if (!_basisVector[i]) {
            r.getElements().push_back(shapelet::ShapeletFunction(0, shapelet::HERMITE, ellipses[i]));
            r.getElements().back().getCoefficients()[0]
                = amplitudes[c] / shapelet::ShapeletFunction::FLUX_FACTOR;
            ++c;
        } else {
            int k = _basisVector[i]->getSize();
            shapelet::MultiShapeletFunction p = _basisVector[i]->makeFunction(
                ellipses[i], amplitudes[ndarray::view(c,c+k)]
            );
            r.getElements().splice(r.getElements().end(), p.getElements());
            c += k;
        }
    }
    return r;
}

Model::Model(BasisVector basisVector, int parameterDim) :
    _parameterDim(parameterDim),
    _amplitudeDim(0),
    _basisVector(basisVector)
{
    for (BasisVector::const_iterator i = _basisVector.begin(); i != _basisVector.end(); ++i) {
        if (!(*i)) {
            ++_amplitudeDim; // null basis indicates a point source
        } else {
            _amplitudeDim += (**i).getSize();
        }
    }
}

// ========== FixedCenterModel ==============================================================================

/*
 * FixedCenterModel uses holds the center fixed and hence has no center parameters.
 * The full parameter vector is ordered
 * [e1[0], e2[0], r[0], e1[1], e2[1], r[1], ...].
 */

namespace {

class FixedCenterModel : public Model {
public:

    FixedCenterModel(BasisVector basisVector, afw::geom::Point2D const & center) :
        Model(basisVector, countParameters(basisVector, 3, 0)), _center(center)
    {}

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const {
        if (prior->getTag() != "single-ellipse") {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                "Cannot adapt prior unless its tag is 'single-ellipse'"
            );
        }
        return prior;
    }

    virtual EllipseVector makeEllipseVector() const {
        return makeEllipseVectorImpl(getBasisVector(), _center);
    }

    virtual void writeEllipses(double const * parameterIter, EllipseIterator ellipseIter) const {
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(parameterIter);
                parameterIter += 3;
            }
        }
    }

    virtual void readEllipses(EllipseConstIterator ellipseIter, double * parameterIter) const {
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().writeParameters(parameterIter);
                parameterIter += 3;
            }
        }
    }

private:
    afw::geom::Point2D _center;
};

} // anonymous

PTR(Model) Model::makeFixedCenter(BasisVector basisVector, afw::geom::Point2D const & center) {
    return boost::make_shared<FixedCenterModel>(basisVector, center);
}

// ========== SingleCenterModel ==============================================================================

/*
 * SingleCenterModel uses one center for all bases.  The full parameter vector is ordered
 * [e1[0], e2[0], r[0], e1[1], e2[1], r[1], ..., x, y].
 */

namespace {

class SingleCenterModel : public Model {
public:

    explicit SingleCenterModel(BasisVector basisVector) :
        Model(basisVector, countParameters(basisVector, 3, 0) + 2)
    {}

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "adaptPrior not implemented for SingleCenterModel"
        );
    }

    virtual EllipseVector makeEllipseVector() const {
        return makeEllipseVectorImpl(getBasisVector(), afw::geom::Point2D());
    }

    virtual void writeEllipses(double const * parameterIter, EllipseIterator ellipseIter) const {
        afw::geom::Point2D center(parameterIter[getParameterDim()-2], parameterIter[getParameterDim()-1]);
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(parameterIter);
                parameterIter += 3;
            }
            ellipseIter->setCenter(center);
        }
    }

    virtual void readEllipses(EllipseConstIterator ellipseIter, double * parameterIter) const {
        // Ellipses have more centers than we need, so we average them.  In most cases, they'll
        // all be the same anyway.
        Eigen::Vector2d p = Eigen::Vector2d::Zero();
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().writeParameters(parameterIter);
                parameterIter += 3;
            }
            p += ellipseIter->getCenter().asEigen();
        }
        p /= getBasisCount();
        parameterIter[0] = p.x(); // note that we've already incremented parameterIter, so these
        parameterIter[1] = p.y(); // are the last two elements, not the first two elements
    }
};

} // anonymous

PTR(Model) Model::makeSingleCenter(BasisVector basisVector) {
    return boost::make_shared<SingleCenterModel>(basisVector);
}

// ========== MultiCenterModel ==============================================================================

/*
 * MultiCenterModel gives each basis its own center parameters.  The full parameter vector is ordered
 * [e1[0], e2[0], r[0], e1[1], e2[1], r[1], ..., x[0], y[0], x[1], y[1], ...].
 */

namespace {

class MultiCenterModel : public Model {
public:

    explicit MultiCenterModel(BasisVector basisVector) :
        Model(basisVector, countParameters(basisVector, 3, 2)),
        _centerParameterOffset(getParameterDim() - getBasisVector().size()*2)
    {}

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "adaptPrior not implemented for MultiCenterModel"
        );
    }

    virtual EllipseVector makeEllipseVector() const {
        return makeEllipseVectorImpl(getBasisVector(), afw::geom::Point2D());
    }

    virtual void writeEllipses(double const * parameterIter, EllipseIterator ellipseIter) const {
        double const * centerIter = parameterIter + _centerParameterOffset;
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(parameterIter);
                parameterIter += 3;
            }
            ellipseIter->setCenter(afw::geom::Point2D(centerIter[0], centerIter[1]));
            centerIter += 2;
        }
    }

    virtual void readEllipses(EllipseConstIterator ellipseIter, double * parameterIter) const {
        double * centerIter = parameterIter + _centerParameterOffset;
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().writeParameters(parameterIter);
                parameterIter += 3;
            }
            centerIter[0] = ellipseIter->getCenter().getX();
            centerIter[1] = ellipseIter->getCenter().getY();
            centerIter += 2;
        }
    }

private:
    int _centerParameterOffset;
};

} // anonymous

PTR(Model) Model::makeMultiCenter(BasisVector basisVector) {
    return boost::make_shared<MultiCenterModel>(basisVector);
}

// ========== MultiModel ====================================================================================

namespace {

static Model::BasisVector concatenateBasisVectors(ModelVector const & components) {
    Model::BasisVector r;
    for (ModelVector::const_iterator i = components.begin(); i != components.end(); ++i) {
        r.insert(r.end(), (**i).getBasisVector().begin(), (**i).getBasisVector().end());
    }
    return r;
}

static int sumParameterDim(ModelVector const & components) {
    int r = 0;
    for (ModelVector::const_iterator i = components.begin(); i != components.end(); ++i) {
        r += (**i).getParameterDim();
    }
    return r;
}

} // anonymous

MultiModel::MultiModel(ModelVector components) :
    Model(concatenateBasisVectors(components), sumParameterDim(components)),
    _components(components)
{}

Model::EllipseVector MultiModel::makeEllipseVector() const {
    EllipseVector r;
    for (ModelVector::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        EllipseVector c = (**i).makeEllipseVector();
        r.insert(r.end(), c.begin(), c.end());
    }
    return r;
}

void MultiModel::writeEllipses(double const * parameterIter, EllipseIterator ellipseIter) const {
    for (ModelVector::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        (**i).writeEllipses(parameterIter, ellipseIter);
        parameterIter += (**i).getParameterDim();
        ellipseIter += (**i).getBasisCount();
    }
}

void MultiModel::readEllipses(EllipseConstIterator ellipseIter, double * parameterIter) const {
    for (ModelVector::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        (**i).readEllipses(ellipseIter, parameterIter);
        parameterIter += (**i).getParameterDim();
        ellipseIter += (**i).getBasisCount();
    }
}

}}} // namespace lsst::meas::multifit
