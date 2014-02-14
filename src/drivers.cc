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

#include "lsst/meas/multifit/drivers.h"

namespace lsst { namespace meas { namespace multifit {

OptimizerFit::OptimizerFit(
    PTR(Model) model, PTR(Prior) prior,
    PTR(afw::coord::Coord) position,
    UnitSystem const & fitSys, UnitSystem const & measSys,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar const,1,1> const & fixed,
    Control const & ctrl
) :
    _model(model), _prior(prior), _position(position), _hasMeasQuantities(false),
    _objectiveValue(std::numeric_limits<Scalar>::quiet_NaN()),
    _flux(std::numeric_limits<Scalar>::quiet_NaN()),
    _fluxSigma(std::numeric_limits<Scalar>::quiet_NaN()),
    _fitSys(fitSys), _measSys(measSys), _fitSysToMeasSys(*position, _fitSys, _measSys),
    _ellipses(_model->makeEllipseVector()),
    _parameters(ndarray::allocate(_model->getNonlinearDim() + model->getAmplitudeDim())),
    _nonlinear(_parameters[ndarray::view(0, _model->getNonlinearDim())]),
    _amplitudes(_parameters[ndarray::view(_model->getNonlinearDim(), _parameters.getSize<0>())]),
    _fixed(ndarray::copy(fixed)),
    _hessian(ndarray::allocate(_parameters.getSize<0>(), _parameters.getSize<0>())),
    _ctrl(ctrl)
{
    LSST_THROW_IF_NE(
        model->getNonlinearDim(), nonlinear.getSize<0>(),
        pex::exceptions::LengthErrorException,
        "Model nonlinear dimension (%d) and nonlinear array size (%d) do not match."
    );
    LSST_THROW_IF_NE(
        model->getAmplitudeDim(), amplitudes.getSize<0>(),
        pex::exceptions::LengthErrorException,
        "Model amplitude dimension (%d) and amplitude array size (%d) do not match."
    );
    LSST_THROW_IF_NE(
        model->getFixedDim(), fixed.getSize<0>(),
        pex::exceptions::LengthErrorException,
        "Model fixed dimension (%d) and fixed array size (%d) do not match."
    );
    _nonlinear.deep() = nonlinear;
    _amplitudes.deep() = amplitudes;
}

void OptimizerFit::run(
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    bool doRecordHistory
) {
    if (doRecordHistory && !_historyRecorder) {
        afw::table::Schema historySchema;
        _historyRecorder.reset(new OptimizerHistoryRecorder(boost::ref(historySchema), _model, true));
        _history = afw::table::BaseCatalog(historySchema);
    }
    PTR(ProjectedLikelihood) likelihood = boost::make_shared<ProjectedLikelihood>(
        _model, _fixed, _fitSys, *_position, exposure, footprint, psf,
        _ctrl.likelihood
    );
    PTR(OptimizerObjective) objective = OptimizerObjective::makeFromLikelihood(likelihood, _prior);
    Optimizer optimizer(objective, _parameters, _ctrl.optimizer);
    _optimizerState = 0;
    if (doRecordHistory) {
        _optimizerState = optimizer.run(*_historyRecorder, _history);
    } else {
        _optimizerState = optimizer.run();
        _history.clear();
    }
    _parameters.deep() = optimizer.getParameters();
    _hessian.deep() = optimizer.getHessian();
    _objectiveValue = optimizer.getObjectiveValue();
    _hasMeasQuantities = false;
}

void OptimizerFit::setModel(PTR(Model) model, bool doForceEllipseConversion) {
    LSST_THROW_IF_NE(
        model->getAmplitudeDim(), _model->getAmplitudeDim(),
        pex::exceptions::LengthErrorException,
        "New model amplitude dimension (%d) does not match old amplitude dimension (%d)"
    );
    if (doForceEllipseConversion
        || model->getNonlinearDim() != _model->getNonlinearDim()
        || model->getFixedDim() != _model->getFixedDim()
    ) {
        if (model->getBasisCount() != _model->getBasisCount()) {
            throw LSST_EXCEPT(
                pex::exceptions::LengthErrorException,
                "The new model must have either the same number of nonlinear and fixed parameters or the"
                " same number of ellipses"
            );
        }
        _model->writeEllipses(_nonlinear.begin(), _fixed.begin(), _ellipses.begin());
        // When we convert parameters via ellipses, we have to make a new EllipseVector for the new Model,
        // because it may assume a different ellipse parametrization than the old one even though it has
        // the same number of ellipses.  When we call std::copy, that invokes assignment operators, which
        // don't change the ellipse type, just its value.
        Model::EllipseVector newEllipses = model->makeEllipseVector();
        _hasMeasQuantities = false;
        std::copy(_ellipses.begin(), _ellipses.end(), newEllipses.begin());
        // Because we allocate the nonlinear and amplitude parameters in one block, if we have to reallocate
        // the former we have to reallocate the latter.
        if (model->getNonlinearDim() != _model->getNonlinearDim()) {
            _parameters = ndarray::allocate(model->getNonlinearDim() + model->getAmplitudeDim());
            _nonlinear = _parameters[ndarray::view(0, model->getNonlinearDim())];
            ndarray::Array<Scalar,1,1> oldAmplitudes = _amplitudes;
            _amplitudes = _parameters[ndarray::view(model->getNonlinearDim(), _parameters.getSize<0>())];
            _amplitudes.deep() = oldAmplitudes;
        }
        if (model->getFixedDim() != _model->getFixedDim()) {
            _fixed = ndarray::allocate(model->getFixedDim());
        }
        model->readEllipses(_ellipses.begin(), _nonlinear.begin(), _fixed.begin());
        _ellipses.swap(newEllipses);
    }
    if (_prior) {
        _prior = model->adaptPrior(_prior);
    }
    _model = model;
}

void OptimizerFit::setPrior(PTR(Prior) prior) {
    if (prior) {
        _prior = _model->adaptPrior(prior);
    } else {
        _prior.reset();
    }
}

void OptimizerFit::setPosition(PTR(afw::coord::Coord) position) {
    _fitSysToMeasSys = LocalUnitTransform(*position, _fitSys, _measSys);
    _position = position;
}

Scalar OptimizerFit::getFlux() const {
    _ensureMeasQuantities();
    return _flux;
}

Scalar OptimizerFit::getFluxSigma() const {
    _ensureMeasQuantities();
    return _fluxSigma;
}

afw::geom::ellipses::Ellipse OptimizerFit::getEllipse(int n) const {
    _ensureMeasQuantities();
    return _ellipses[n];
}

void OptimizerFit::_ensureMeasQuantities() const {
    if (!_hasMeasQuantities) {
        _model->writeEllipses(_nonlinear.begin(), _fixed.begin(), _ellipses.begin());
        for (Model::EllipseVector::iterator iter = _ellipses.begin(); iter != _ellipses.end(); ++iter) {
            iter->transform(_fitSysToMeasSys.geometric).inPlace();
        }
        _flux = _amplitudes.asEigen().sum() * _fitSysToMeasSys.flux;
        _fluxSigma = std::numeric_limits<Scalar>::quiet_NaN(); // tricky; will do later
        _hasMeasQuantities = true;
    }
}

}}} // namespace lsst::meas::multifit
