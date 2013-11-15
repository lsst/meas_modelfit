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
#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

static afw::geom::ellipses::SeparableConformalShearLogTraceRadius extendedSourceEllipseCore;
static afw::geom::ellipses::Quadrupole pointSourceEllipseCore(0.0, 0.0, 0.0);

int countNonlinears(Model::BasisVector const & basisVector, int ellipseFactor, int pointFactor) {
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

Model::EllipseVector makeEllipseVectorImpl(Model::BasisVector const & basisVector) {
    Model::EllipseVector r;
    r.reserve(basisVector.size());
    for (int i = 0, n = basisVector.size(); i < n; ++i) {
        if (basisVector[i]) {
            r.push_back(afw::geom::ellipses::Ellipse(extendedSourceEllipseCore));
        } else {
            r.push_back(afw::geom::ellipses::Ellipse(pointSourceEllipseCore));
        }
    }
    return r;
}

} // anonymous

// ========== Model =========================================================================================

shapelet::MultiShapeletFunction Model::makeShapeletFunction(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar const,1,1> const & fixed
) const {
    EllipseVector ellipses = makeEllipseVector();
    writeEllipses(nonlinear.begin(), fixed.begin(), ellipses.begin());
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

Model::Model(
    BasisVector basisVector,
    NameVector nonlinearNames,
    NameVector amplitudeNames,
    NameVector fixedNames
) :
    _nonlinearNames(nonlinearNames),
    _amplitudeNames(amplitudeNames),
    _fixedNames(fixedNames),
    _basisVector(basisVector)
{
    int amplitudeDim = 0;
    for (BasisVector::const_iterator i = _basisVector.begin(); i != _basisVector.end(); ++i) {
        if (!(*i)) {
            ++amplitudeDim; // null basis indicates a point source
        } else {
            amplitudeDim += (**i).getSize();
        }
    }
    LSST_ASSERT_EQUAL(
        amplitudeDim, getAmplitudeDim(),
        "Number of amplitudes in basis vectors (%d) does not match number of amplitude names (%d)",
        pex::exceptions::LengthErrorException
    );
}

Model::EllipseVector Model::writeEllipses(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & fixed
) const {
    LSST_ASSERT_EQUAL(
        nonlinear.getSize<0>(), getNonlinearDim(),
        "Size of nonlinear array (%d) does not match dimension of model (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        fixed.getSize<0>(), getFixedDim(),
        "Size of fixed array (%d) does not match dimension of model (%d)",
        pex::exceptions::LengthErrorException
    );
    EllipseVector r = makeEllipseVector();
    writeEllipses(nonlinear.begin(), fixed.begin(), r.begin());
    return r;
}

void Model::readEllipses(
    EllipseVector const & ellipses,
    ndarray::Array<Scalar,1,1> const & nonlinear,
    ndarray::Array<Scalar,1,1> const & fixed
) const {
    LSST_ASSERT_EQUAL(
        nonlinear.getSize<0>(), getNonlinearDim(),
        "Size of nonlinear array (%d) does not match dimension of model (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        fixed.getSize<0>(), getFixedDim(),
        "Size of fixed array (%d) does not match dimension of model (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        int(ellipses.size()), getBasisCount(),
        "Size of ellipse vector (%d) does not match basis count (%d)",
        pex::exceptions::LengthErrorException
    );
    readEllipses(ellipses.begin(), nonlinear.begin(), fixed.begin());
}

// ========== FixedCenterModel ==============================================================================

/*
 * FixedCenterModel uses holds the center fixed and hence has no free center parameters.
 * The nonlinear vector is ordered [e1[0], e2[0], r[0], e1[1], e2[1], r[1], ...], while
 * the fixed parameter vector is simply [x,y].
 */

namespace {

class FixedCenterModel : public Model {
public:

    static NameVector getFixedNames() {
        NameVector r(2);
        r[0] = "x";
        r[1] = "y";
        return r;
    }

    FixedCenterModel(BasisVector basisVector, NameVector nonlinearNames, NameVector amplitudeNames) :
        Model(basisVector, nonlinearNames, amplitudeNames, getFixedNames())
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
        return makeEllipseVectorImpl(getBasisVector());
    }

    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const {
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(nonlinearIter);
                ellipseIter->getCenter().setX(fixedIter[0]);
                ellipseIter->getCenter().setY(fixedIter[1]);
                nonlinearIter += 3;
            }
        }
    }

    virtual void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const {
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().writeParameters(nonlinearIter);
                fixedIter[0] = ellipseIter->getCenter().getX();
                fixedIter[1] = ellipseIter->getCenter().getY();
                nonlinearIter += 3;
            }
        }
    }

};

} // anonymous

PTR(Model) Model::makeFixedCenter(BasisVector basisVector, NameVector const & prefixes) {
    LSST_ASSERT_EQUAL(
        basisVector.size(), prefixes.size(),
        "Size of basis vector (%d) does not match number of prefixes (%d)",
        pex::exceptions::LengthErrorException
    );
    NameVector nonlinearNames;
    NameVector amplitudeNames;
    for (std::size_t i = 0, n = basisVector.size(); i < n; ++i) {
        if (basisVector[i]) {
            nonlinearNames.push_back(prefixes[i] + "eta1");
            nonlinearNames.push_back(prefixes[i] + "eta2");
            nonlinearNames.push_back(prefixes[i] + "logR");
        }
        for (std::size_t j = 0, m = basisVector[i]->getSize(); j < m; ++j) {
            amplitudeNames.push_back((boost::format("%salpha%d") % prefixes[i] % j).str());
        }
    }
    return boost::make_shared<FixedCenterModel>(basisVector, nonlinearNames, amplitudeNames);
}

PTR(Model) Model::makeFixedCenter(PTR(shapelet::MultiShapeletBasis) basis) {
    NameVector nonlinearNames;
    NameVector amplitudeNames;
    if (basis) {
        nonlinearNames.push_back("eta1");
        nonlinearNames.push_back("eta2");
        nonlinearNames.push_back("logR");
    }
    for (std::size_t j = 0, m = basis->getSize(); j < m; ++j) {
        amplitudeNames.push_back((boost::format("alpha%d") % j).str());
    }
    BasisVector basisVector(1, basis);
    return boost::make_shared<FixedCenterModel>(basisVector, nonlinearNames, amplitudeNames);
}

// ========== SingleCenterModel ==============================================================================

/*
 * SingleCenterModel uses one center for all bases.  The nonlinear vector is ordered
 * [e1[0], e2[0], r[0], e1[1], e2[1], r[1], ..., x, y], while the fixed parameter is empty.
 */

namespace {

class SingleCenterModel : public Model {
public:

    SingleCenterModel(BasisVector basisVector, NameVector nonlinearNames, NameVector amplitudeNames) :
        Model(basisVector, nonlinearNames, amplitudeNames, NameVector())
    {}

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "adaptPrior not implemented for SingleCenterModel"
        );
    }

    virtual EllipseVector makeEllipseVector() const {
        return makeEllipseVectorImpl(getBasisVector());
    }

    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const {
        afw::geom::Point2D center(nonlinearIter[getNonlinearDim()-2], nonlinearIter[getNonlinearDim()-1]);
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(nonlinearIter);
                nonlinearIter += 3;
            }
            ellipseIter->setCenter(center);
        }
    }

    virtual void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const {
        // Ellipses have more centers than we need, so we average them.  In most cases, they'll
        // all be the same anyway.
        Eigen::Vector2d p = Eigen::Vector2d::Zero();
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().writeParameters(nonlinearIter);
                nonlinearIter += 3;
            }
            p += ellipseIter->getCenter().asEigen();
        }
        p /= getBasisCount();
        nonlinearIter[0] = p.x(); // note that we've already incremented nonlinearIter, so these
        nonlinearIter[1] = p.y(); // are the last two elements, not the first two elements
    }
};

} // anonymous

PTR(Model) Model::makeSingleCenter(BasisVector basisVector, NameVector const & prefixes) {
    LSST_ASSERT_EQUAL(
        basisVector.size(), prefixes.size(),
        "Size of basis vector (%d) does not match number of prefixes (%d)",
        pex::exceptions::LengthErrorException
    );
    NameVector nonlinearNames;
    NameVector amplitudeNames;
    for (std::size_t i = 0, n = basisVector.size(); i < n; ++i) {
        if (basisVector[i]) {
            nonlinearNames.push_back(prefixes[i] + "eta1");
            nonlinearNames.push_back(prefixes[i] + "eta2");
            nonlinearNames.push_back(prefixes[i] + "logR");
        }
        for (std::size_t j = 0, m = basisVector[i]->getSize(); j < m; ++j) {
            amplitudeNames.push_back((boost::format("%salpha%d") % prefixes[i] % j).str());
        }
    }
    nonlinearNames.push_back("x");
    nonlinearNames.push_back("y");
    return boost::make_shared<SingleCenterModel>(basisVector, nonlinearNames, amplitudeNames);
}

PTR(Model) Model::makeSingleCenter(PTR(shapelet::MultiShapeletBasis) basis) {
    NameVector nonlinearNames;
    NameVector amplitudeNames;
    if (basis) {
        nonlinearNames.push_back("eta1");
        nonlinearNames.push_back("eta2");
        nonlinearNames.push_back("logR");
    }
    for (std::size_t j = 0, m = basis->getSize(); j < m; ++j) {
        amplitudeNames.push_back((boost::format("alpha%d") % j).str());
    }
    nonlinearNames.push_back("x");
    nonlinearNames.push_back("y");
    BasisVector basisVector(1, basis);
    return boost::make_shared<FixedCenterModel>(basisVector, nonlinearNames, amplitudeNames);
}

// ========== MultiCenterModel ==============================================================================

/*
 * MultiCenterModel gives each basis its own center parameters.  The full nonlinear vector is ordered
 * [e1[0], e2[0], r[0], e1[1], e2[1], r[1], ..., x[0], y[0], x[1], y[1], ...].
 */

namespace {

class MultiCenterModel : public Model {
public:

    explicit MultiCenterModel(BasisVector basisVector, NameVector nonlinearNames, NameVector amplitudeNames) :
        Model(basisVector, nonlinearNames, amplitudeNames, NameVector()),
        _centerParameterOffset(getNonlinearDim() - getBasisVector().size()*2)
    {}

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "adaptPrior not implemented for MultiCenterModel"
        );
    }

    virtual EllipseVector makeEllipseVector() const {
        return makeEllipseVectorImpl(getBasisVector());
    }

    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const {
        Scalar const * centerIter = nonlinearIter + _centerParameterOffset;
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(nonlinearIter);
                nonlinearIter += 3;
            }
            ellipseIter->setCenter(afw::geom::Point2D(centerIter[0], centerIter[1]));
            centerIter += 2;
        }
    }

    virtual void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const {
        Scalar * centerIter = nonlinearIter + _centerParameterOffset;
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().writeParameters(nonlinearIter);
                nonlinearIter += 3;
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

PTR(Model) Model::makeMultiCenter(BasisVector basisVector, NameVector const & prefixes) {
    LSST_ASSERT_EQUAL(
        basisVector.size(), prefixes.size(),
        "Size of basis vector (%d) does not match number of prefixes (%d)",
        pex::exceptions::LengthErrorException
    );
    NameVector nonlinearNames;
    NameVector amplitudeNames;
    for (std::size_t i = 0, n = basisVector.size(); i < n; ++i) {
        if (basisVector[i]) {
            nonlinearNames.push_back(prefixes[i] + "eta1");
            nonlinearNames.push_back(prefixes[i] + "eta2");
            nonlinearNames.push_back(prefixes[i] + "logR");
        }
        for (std::size_t j = 0, m = basisVector[i]->getSize(); j < m; ++j) {
            amplitudeNames.push_back((boost::format("%salpha%d") % prefixes[i] % j).str());
        }
    }
    for (std::size_t i = 0, n = basisVector.size(); i < n; ++i) {
        nonlinearNames.push_back(prefixes[i] + "x");
        nonlinearNames.push_back(prefixes[i] + "y");
    }
    return boost::make_shared<MultiCenterModel>(basisVector, nonlinearNames, amplitudeNames);
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

typedef Model::NameVector const & (Model::*ModelNameGetter)() const;

static Model::NameVector concatenateNameVectors(
    ModelVector const & components, Model::NameVector const & prefixes, ModelNameGetter getter
) {
    LSST_ASSERT_EQUAL(
        components.size(), prefixes.size(),
        "Number of model components (%d) does not match number of prefixes (%d)",
        pex::exceptions::LengthErrorException
    );
    Model::NameVector r;
    for (ModelVector::const_iterator i = components.begin(); i != components.end(); ++i) {
        Model::NameVector const & componentNames = ((**i).*getter)();
        for (std::size_t j = 0, n = componentNames.size(); j < n; ++j) {
            r.push_back(prefixes[j] + componentNames[j]);
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

PTR(Prior) MultiModel::adaptPrior(PTR(Prior) prior) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
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

}}} // namespace lsst::meas::multifit
