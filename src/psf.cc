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
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/psf.h"

namespace lsst { namespace meas { namespace multifit {

PTR(Model) makeMultiShapeletPsfModel(std::vector<int> const & orders) {
    ModelVector components;
    Model::BasisVector basisVector(1);
    double radius = -orders.size() * 0.5;
    int fullBasisSize = 0;
    for (std::size_t i = 0; i != orders.size(); ++i) {
        fullBasisSize += shapelet::computeSize(orders[i]);
    }
    int offset = 0;
    PTR(shapelet::MultiShapeletBasis) & basis = basisVector.front();
    for (std::size_t i = 0; i != orders.size(); ++i) {
        int componentSize = shapelet::computeSize(orders[i]);
        ndarray::Array<double,2,2> matrix = ndarray::allocate(componentSize, fullBasisSize);
        matrix.deep() = 0.0;
        matrix.asEigen().block(0, componentSize, offset, componentSize).setIdentity();
        basis = boost::make_shared<shapelet::MultiShapeletBasis>(fullBasisSize);
        basis->addComponent(radius, orders[i], matrix);
        components.push_back(Model::makeMultiCenter(basisVector));
    }
    return boost::make_shared<MultiModel>(components);
}

class PsfApproximationLikelihood::Impl {
public:

    explicit Impl(
        ndarray::Array<Pixel const,1,1> const & x,
        ndarray::Array<Pixel const,1,1> const & y,
        Model::EllipseVector const & ellipses,
        Model::BasisVector const & basisVector
    ) : _ellipses(ellipses),
        _builder(x, y)
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
    }

private:
    Model::EllipseVector _ellipses;
    shapelet::ModelBuilder<Pixel> _builder;
    ndarray::Array<Pixel,2,-2> _workspace;
};

PsfApproximationLikelihood::PsfApproximationLikelihood(
    ndarray::Array<Pixel const,2,2> const & image,
    afw::geom::Extent2I const & xy0,
    PTR(Model) model,
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
    _impl.reset(new Impl(x, y, model->makeEllipseVector(), model->getBasisVector()));
    _data = ndarray::copy(ndarray::flatten<1>(image));
}

void PsfApproximationLikelihood::computeModelMatrix(
    ndarray::Array<Pixel,2,-1> const & modelMatrix,
    ndarray::Array<Scalar const,1,1> const & nonlinear
) const {
    return _impl->computeModelMatrix(modelMatrix, nonlinear, _fixed, *getModel());
}

PsfApproximationLikelihood::~PsfApproximationLikelihood() {}

}}} // namespace lsst::meas::multifit
