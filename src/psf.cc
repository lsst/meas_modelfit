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
#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/psf.h"

namespace lsst { namespace meas { namespace multifit {


PsfFitter::PsfFitter(PsfFitterControl const & ctrl) : _ctrl(ctrl) {}

PTR(Model) PsfFitter::getInitialModel() const { return _initialModel; }

PTR(Model) PsfFitter::getFinalModel() const { return _finalModel; }

PTR(Prior) PsfFitter::getInitialPrior() const { return _initialPrior; }

PTR(Prior) PsfFitter::getFinalPrior() const { return _finalPrior; }

shapelet::MultiShapeletFunction PsfFitter::fitInitial(
    afw::image::Image<Pixel> const & image,
    Scalar noiseSigma
) {
    shapelet::MultiShapeletFunction result;
    return result;
}

shapelet::MultiShapeletFunction PsfFitter::fitFinal(
    afw::image::Image<Pixel> const & image,
    Scalar noiseSigma,
    shapelet::MultiShapeletFunction const & initialResult
) {
    shapelet::MultiShapeletFunction result;
    return result;
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
        _builders(),
        _sigma(sigma)
    {
        FactoryVector factories;
        factories.reserve(basisVector.size());
        _builders.reserve(basisVector.size());
        int workspaceSize = 0;
        for (Model::BasisVector::const_iterator i = basisVector.begin(); i != basisVector.end(); ++i) {
            factories.push_back(shapelet::MatrixBuilderFactory<Pixel>(x, y, **i));
            workspaceSize = std::max(workspaceSize, factories.back().computeWorkspace());
        }
        shapelet::MatrixBuilderWorkspace<Pixel> workspace(workspaceSize);
        for (FactoryVector::const_iterator i = factories.begin(); i != factories.end(); ++i) {
            shapelet::MatrixBuilderWorkspace<Pixel> wsCopy(workspace); // share workspace between builders
            _builders.push_back((*i)(wsCopy));
        }
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
            _builders[i](modelMatrix, _ellipses[i]);
        }
        modelMatrix.asEigen() /= _sigma;
    }

private:
    typedef std::vector< shapelet::MatrixBuilder<Pixel> > BuilderVector;
    typedef std::vector< shapelet::MatrixBuilderFactory<Pixel> > FactoryVector;

    Model::EllipseVector _ellipses;
    BuilderVector _builders;
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
