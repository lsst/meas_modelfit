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

#include "ndarray/eigen.h"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/modelfit/Model.h"
#include "lsst/meas/modelfit/Prior.h"
#include "lsst/meas/modelfit/UnitSystem.h"

namespace lsst { namespace meas { namespace modelfit {

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

    std::shared_ptr<Prior> adaptPrior(std::shared_ptr<Prior> prior) const override {
        if (prior->getTag() != "single-ellipse") {
            throw LSST_EXCEPT(
                pex::exceptions::LogicError,
                "Cannot adapt prior unless its tag is 'single-ellipse'"
            );
        }
        return prior;
    }

    EllipseVector makeEllipseVector() const override {
        return makeEllipseVectorImpl(getBasisVector());
    }

    void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const override {
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(nonlinearIter);
                ellipseIter->getCenter().setX(fixedIter[0]);
                ellipseIter->getCenter().setY(fixedIter[1]);
                nonlinearIter += 3;
            }
        }
    }

    void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const override {
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

    std::shared_ptr<Prior> adaptPrior(std::shared_ptr<Prior> prior) const override {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "adaptPrior not implemented for SingleCenterModel"
        );
    }

    EllipseVector makeEllipseVector() const override {
        return makeEllipseVectorImpl(getBasisVector());
    }

    void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const override {
        geom::Point2D center(nonlinearIter[getNonlinearDim()-2], nonlinearIter[getNonlinearDim()-1]);
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(nonlinearIter);
                nonlinearIter += 3;
            }
            ellipseIter->setCenter(center);
        }
    }

    void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const override {
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

    std::shared_ptr<Prior> adaptPrior(std::shared_ptr<Prior> prior) const override {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "adaptPrior not implemented for MultiCenterModel"
        );
    }

    EllipseVector makeEllipseVector() const override {
        return makeEllipseVectorImpl(getBasisVector());
    }

    void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter,
        EllipseIterator ellipseIter
    ) const override {
        Scalar const * centerIter = nonlinearIter + _centerParameterOffset;
        for (int i = 0; i < getBasisCount(); ++i, ++ellipseIter) {
            if (getBasisVector()[i]) {
                ellipseIter->getCore().readParameters(nonlinearIter);
                nonlinearIter += 3;
            }
            ellipseIter->setCenter(geom::Point2D(centerIter[0], centerIter[1]));
            centerIter += 2;
        }
    }

    void readEllipses(
        EllipseConstIterator ellipseIter,
        Scalar * nonlinearIter, Scalar * fixedIter
    ) const override {
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
            r.getComponents().push_back(shapelet::ShapeletFunction(0, shapelet::HERMITE, ellipses[i]));
            r.getComponents().back().getCoefficients()[0]
                = amplitudes[c] / shapelet::ShapeletFunction::FLUX_FACTOR;
            ++c;
        } else {
            int k = _basisVector[i]->getSize();
            shapelet::MultiShapeletFunction p = _basisVector[i]->makeFunction(
                ellipses[i], amplitudes[ndarray::view(c,c+k)]
            );
            r.getComponents().insert(r.getComponents().end(), p.getComponents().begin(), p.getComponents().end());
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
    LSST_THROW_IF_NE(
        amplitudeDim, getAmplitudeDim(),
        pex::exceptions::LengthError,
        "Number of amplitudes in basis vectors (%d) does not match number of amplitude names (%d)"
    );
}

Model::EllipseVector Model::writeEllipses(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & fixed
) const {
    LSST_THROW_IF_NE(
        nonlinear.getSize<0>(), getNonlinearDim(),
        pex::exceptions::LengthError,
        "Size of nonlinear array (%d) does not match dimension of model (%d)"
    );
    LSST_THROW_IF_NE(
        fixed.getSize<0>(), getFixedDim(),
        pex::exceptions::LengthError,
        "Size of fixed array (%d) does not match dimension of model (%d)"
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
    LSST_THROW_IF_NE(
        nonlinear.getSize<0>(), getNonlinearDim(),
        pex::exceptions::LengthError,
        "Size of nonlinear array (%d) does not match dimension of model (%d)"
    );
    LSST_THROW_IF_NE(
        fixed.getSize<0>(), getFixedDim(),
        pex::exceptions::LengthError,
        "Size of fixed array (%d) does not match dimension of model (%d)"
    );
    LSST_THROW_IF_NE(
        int(ellipses.size()), getBasisCount(),
        pex::exceptions::LengthError,
        "Size of ellipse vector (%d) does not match basis count (%d)"
    );
    readEllipses(ellipses.begin(), nonlinear.begin(), fixed.begin());
}

void Model::transformParameters(
    LocalUnitTransform const & transform,
    ndarray::Array<Scalar,1,1> const & nonlinear,
    ndarray::Array<Scalar,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & fixed
) const {
    EllipseVector ellipses = writeEllipses(nonlinear, fixed);
    for (EllipseVector::iterator i = ellipses.begin(); i != ellipses.end(); ++i) {
        i->transform(transform.geometric).inPlace();
    }
    readEllipses(ellipses.begin(), nonlinear.begin(), fixed.begin());
    ndarray::asEigenMatrix(amplitudes) *= transform.flux;
}


std::shared_ptr<Model> Model::make(BasisVector basisVector, NameVector const & prefixes, CenterEnum center) {
    LSST_THROW_IF_NE(
        basisVector.size(), prefixes.size(),
        pex::exceptions::LengthError,
        "Size of basis vector (%d) does not match number of prefixes (%d)"
    );
    NameVector nonlinearNames;
    NameVector amplitudeNames;
    if (center == FIXED_CENTER || center == SINGLE_CENTER) {
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
        if (center == FIXED_CENTER) {
            return std::make_shared<FixedCenterModel>(basisVector, nonlinearNames, amplitudeNames);
        } else {
            nonlinearNames.push_back("x");
            nonlinearNames.push_back("y");
            return std::make_shared<SingleCenterModel>(basisVector, nonlinearNames, amplitudeNames);
        }
    } else if (center == MULTI_CENTER) {
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
        return std::make_shared<MultiCenterModel>(basisVector, nonlinearNames, amplitudeNames);
    }
    throw LSST_EXCEPT(
        pex::exceptions::LogicError,
        "Unexpected value for CenterEnum"
    );
}

std::shared_ptr<Model> Model::make(std::shared_ptr<shapelet::MultiShapeletBasis> basis, CenterEnum center) {
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
    if (center == FIXED_CENTER) {
        BasisVector basisVector(1, basis);
        return std::make_shared<FixedCenterModel>(basisVector, nonlinearNames, amplitudeNames);
    } else if (center == SINGLE_CENTER || center == MULTI_CENTER) {
        nonlinearNames.push_back("x");
        nonlinearNames.push_back("y");
        BasisVector basisVector(1, basis);
        return std::make_shared<SingleCenterModel>(basisVector, nonlinearNames, amplitudeNames);
    }
    throw LSST_EXCEPT(
        pex::exceptions::LogicError,
        "Unexpected value for CenterEnum"
    );
}

std::shared_ptr<Model> Model::makeGaussian(CenterEnum center, double radius) {
    std::shared_ptr<shapelet::MultiShapeletBasis> basis = std::make_shared<shapelet::MultiShapeletBasis>(1);
    ndarray::Array<double,2,2> matrix = ndarray::allocate(1, 1);
    matrix[0][0] = 1.0;
    basis->addComponent(radius, 0, matrix);
    basis->normalize();
    return make(basis, center);
}

}}} // namespace lsst::meas::modelfit
