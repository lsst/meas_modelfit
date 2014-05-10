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
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/shapelet/GaussHermiteEvaluator.h" // for PackedIndex
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/psf.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

typedef std::vector< std::pair<std::string,PsfFitterComponentControl> > ComponentVector;
typedef ComponentVector::const_iterator ComponentIterator;

static afw::geom::ellipses::SeparableConformalShearLogTraceRadius psfFitterEllipseCore;

class PsfFitterModel : public Model {
public:

    PsfFitterModel(
        BasisVector basisVector,
        NameVector nonlinearNames,
        NameVector amplitudeNames,
        NameVector fixedNames,
        ComponentVector components
    ) : Model(basisVector, nonlinearNames, amplitudeNames, fixedNames), _components(components) {}

    virtual PTR(Prior) adaptPrior(PTR(Prior) prior) const {
        return prior; // TODO
    }

    virtual EllipseVector makeEllipseVector() const {
        return EllipseVector(getBasisCount(), afw::geom::ellipses::Ellipse(psfFitterEllipseCore));
    }

    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter, EllipseIterator ellipseIter
    ) const {
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i, ++ellipseIter) {
            afw::geom::ellipses::Ellipse::ParameterVector v;
            // we always store fiducial values in the fixed parameter array; if a parameter is allowed
            // to vary, we store the offset from the fiducial value in the nonlinear parameter array.
            for (int n = 0; n < 5; ++n, ++fixedIter) {
                v[n] = *fixedIter;
            }
            if (i->second.ellipticityPriorSigma > 0.0) {
                v[0] += *(nonlinearIter++);
                v[1] += *(nonlinearIter++);
            }
            if (i->second.radiusPriorSigma > 0.0) {
                v[2] += *(nonlinearIter++);
            }
            if (i->second.positionPriorSigma > 0.0) {
                v[3] += *(nonlinearIter++);
                v[4] += *(nonlinearIter++);
            }
            ellipseIter->setParameterVector(v);
        }
    }

    virtual void readEllipses(
        EllipseConstIterator ellipseIter, Scalar * nonlinearIter, Scalar * fixedIter
    ) const {
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i, ++ellipseIter) {
            afw::geom::ellipses::Ellipse::ParameterVector v = ellipseIter->getParameterVector();
            // We always store fiducial values in the fixed parameter array; if a parameter is allowed
            // to vary, we store the offset from the fiducial value in the nonlinear parameter array.
            // When reading ellipses, we only update the fixed array, and set the offsets in the
            // nonlinear array to zero.
            // I'm a little concerned that this means we can't round-trip parameter vectors through
            // ellipses with this model, but I can't think of a concrete case where that would cause
            // problems.
            for (int n = 0; n < 5; ++n, ++fixedIter) {
                *fixedIter = v[n];
            }
            if (i->second.ellipticityPriorSigma > 0.0) {
                *(nonlinearIter++) = 0;
                *(nonlinearIter++) = 0;
            }
            if (i->second.radiusPriorSigma > 0.0) {
                *(nonlinearIter++) = 0;
            }
            if (i->second.positionPriorSigma > 0.0) {
                *(nonlinearIter++) = 0;
                *(nonlinearIter++) = 0;
            }
        }
    }

private:
    ComponentVector _components;
};

ComponentVector vectorizeComponents(PsfFitterControl const & ctrl) {
    ComponentVector components;
    if (ctrl.inner.order >= 0) {
        components.push_back(std::make_pair("inner", ctrl.inner));
    }
    if (ctrl.primary.order >= 0) {
        components.push_back(std::make_pair("primary", ctrl.primary));
    }
    if (ctrl.wings.order >= 0) {
        components.push_back(std::make_pair("wings", ctrl.wings));
    }
    if (ctrl.outer.order >= 0) {
        components.push_back(std::make_pair("outer", ctrl.outer));
    }
    return components;
}

} // anonymous

PsfFitter::PsfFitter(PsfFitterControl const & ctrl) :
    _ctrl(ctrl)
{
    ComponentVector components = vectorizeComponents(_ctrl);
    Model::BasisVector basisVector;
    Model::NameVector nonlinearNames;
    Model::NameVector amplitudeNames;
    Model::NameVector fixedNames;

    for (ComponentIterator i = components.begin(); i != components.end(); ++i) {

        // Construct a MultiShapeletBasis with a single shapelet basis with this component
        int dim = shapelet::computeSize(i->second.order);
        PTR(shapelet::MultiShapeletBasis) basis = boost::make_shared<shapelet::MultiShapeletBasis>(dim);
        ndarray::Array<double,2,2> matrix = ndarray::allocate(dim, dim);
        matrix.asEigen().setIdentity();
        basis->addComponent(1.0, i->second.order, matrix);
        basisVector.push_back(basis);

        // Append to the name vectors for this component; all of this has to be consistent with the
        // iteration that happens in PsfFitterModel and PsfFitterPrior
        for (shapelet::PackedIndex s; s.getOrder() <= i->second.order; ++s) {
            amplitudeNames.push_back(
                (boost::format("%s.alpha[%d,%d]") % i->first % s.getX() % s.getY()).str()
            );
        }
        fixedNames.push_back(i->first + ".fiducial.eta1");
        fixedNames.push_back(i->first + ".fiducial.eta2");
        fixedNames.push_back(i->first + ".fiducial.logR");
        fixedNames.push_back(i->first + ".fiducial.x");
        fixedNames.push_back(i->first + ".fiducial.y");
        if (i->second.ellipticityPriorSigma > 0.0) {
            nonlinearNames.push_back(i->first + ".eta1");
            nonlinearNames.push_back(i->first + ".eta2");
        }
        if (i->second.radiusPriorSigma > 0.0) {
            nonlinearNames.push_back(i->first + ".logR");
        }
        if (i->second.positionPriorSigma > 0.0) {
            nonlinearNames.push_back(i->first + ".x");
            nonlinearNames.push_back(i->first + ".y");
        }

    }

    _model = boost::make_shared<PsfFitterModel>(
        basisVector, nonlinearNames, amplitudeNames, fixedNames, components
    );
}

shapelet::MultiShapeletFunction PsfFitter::adapt(
    shapelet::MultiShapeletFunction const & previousFit,
    PTR(Model) previousModel
) const {
    
}


shapelet::MultiShapeletFunction PsfFitter::apply(
    afw::image::Image<Pixel> const & image,
    Scalar noiseSigma,
    afw::geom::ellipses::Quadrupole const & moments
) const {

}

shapelet::MultiShapeletFunction PsfFitter::apply(
    afw::image::Image<Pixel> const & image,
    Scalar noiseSigma,
    shapelet::MultiShapeletFunction const & initial
) const {

}


class MultiShapeletPsfLikelihood::Impl {
public:

    explicit Impl(
        ndarray::Array<Pixel const,1,1> const & x,
        ndarray::Array<Pixel const,1,1> const & y,
        Model::EllipseVector const & ellipses,
        Model::BasisVector const & basisVector,
        Scalar sigma
    ) : _ellipses(ellipses),
        _builder(x, y),
        _sigma(sigma)
    {
        int maxOrder = 0;
        for (Model::BasisVector::const_iterator i = basisVector.begin(); i != basisVector.end(); ++i) {
            for (shapelet::MultiShapeletBasis::Iterator j = (**i).begin(); j != (**i).end(); ++j) {
                maxOrder = std::max(j->getOrder(), maxOrder);
            }
        }
        ndarray::Array<Pixel,2,2> workspaceT
            = ndarray::allocate(shapelet::computeSize(maxOrder), x.getSize<0>());
        _workspace = workspaceT.transpose();
    }

    void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & fixed,
        Model const & model
    ) {
        model.writeEllipses(nonlinear.begin(), fixed.begin(), _ellipses.begin());
        modelMatrix.deep() = 0.0;
        Model::BasisVector const & basisVector = model.getBasisVector();
        for (std::size_t i = 0; i < basisVector.size(); ++i) {
            _builder.update(_ellipses[i]);
            for (shapelet::MultiShapeletBasis::Iterator j = basisVector[i]->begin();
                 j != basisVector[i]->end(); ++j
            ) {
                _workspace.deep() = 0.0;
                _builder.addModelMatrix(j->getOrder(), _workspace);
                modelMatrix.asEigen() += _workspace.asEigen() * j->getMatrix().asEigen().cast<Pixel>();
            }
        }
        modelMatrix.asEigen() /= _sigma;
    }

private:
    Model::EllipseVector _ellipses;
    shapelet::ModelBuilder<Pixel> _builder;
    ndarray::Array<Pixel,2,-2> _workspace;
    Scalar _sigma;
};

MultiShapeletPsfLikelihood::MultiShapeletPsfLikelihood(
    ndarray::Array<Pixel const,2,2> const & image,
    afw::geom::Point2I const & xy0,
    PTR(Model) model,
    Scalar sigma,
    ndarray::Array<Scalar const,1,1> const & fixed
) :
    Likelihood(model, fixed)
{
    int nx = image.getSize<1>();
    int ny = image.getSize<0>();
    int nTot = nx*ny;
    ndarray::Array<Pixel,1,1> x = ndarray::allocate(nTot);
    ndarray::Array<Pixel,1,1> y = ndarray::allocate(nTot);
    int j = 0;
    for (int iy = xy0.getY(), yEnd = xy0.getY() + ny; iy < yEnd; ++iy) {
        for (int ix = xy0.getX(), xEnd = xy0.getX() + nx; ix < xEnd; ++ix, ++j) {
            x[j] = ix;
            y[j] = iy;
        }
    }
    _impl.reset(new Impl(x, y, model->makeEllipseVector(), model->getBasisVector(), sigma));
    _data = ndarray::copy(ndarray::flatten<1>(image));
    _data.deep() /= sigma;
    _weights = ndarray::allocate(_data.getShape());
    _weights.deep() = 1.0;
}

void MultiShapeletPsfLikelihood::computeModelMatrix(
    ndarray::Array<Pixel,2,-1> const & modelMatrix,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    bool doApplyWeights
) const {
    return _impl->computeModelMatrix(modelMatrix, nonlinear, _fixed, *getModel());
}

MultiShapeletPsfLikelihood::~MultiShapeletPsfLikelihood() {}

}}} // namespace lsst::meas::multifit
