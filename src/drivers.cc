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

OptimizerFitter::OptimizerFitter(PTR(Model) model, PTR(Prior) prior, Control const & ctrl) :
    _model(model),
    _prior(prior),
    _parameters(ndarray::allocate(_model->getNonlinearDim() + _model->getAmplitudeDim())),
    _ctrl(ctrl)
{
    if (_ctrl.doRecordHistory) {
        _historyRecorder.reset(new OptimizerHistoryRecorder(_historySchema, _model, true));
    }
}

ndarray::Array<Scalar const,1,1> OptimizerFitter::getDefaultInitialAmplitudes() const {
    ndarray::Array<Scalar,1,1> r(ndarray::allocate(_model->getAmplitudeDim()));
    r.deep() = 0.0;
    r[0] = 1.0;
    return r;
}

OptimizerFitter::Result OptimizerFitter::apply(
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    afw::coord::Coord const & position,
    UnitSystem const & measSys,
    UnitSystem const & fitSys,
    ndarray::Array<Scalar const,1,1> const & initialNonlinear,
    ndarray::Array<Scalar const,1,1> const & initialAmplitudes,
    ndarray::Array<Scalar const,1,1> const & fixed,
    shapelet::MultiShapeletFunction const & psf
) const {
    int const nonlinearDim = _model->getNonlinearDim();
    int const amplitudeDim = _model->getAmplitudeDim();
    int const parameterDim = nonlinearDim + amplitudeDim;
    Result result(_model, position, measSys, fitSys);
    result._model = _model;
    result._fixed = fixed;
    PTR(ProjectedLikelihood) likelihood = boost::make_shared<ProjectedLikelihood>(
        _model, fixed, fitSys, position, exposure, footprint, psf,
        _ctrl.likelihood
    );
    PTR(OptimizerObjective) objective = OptimizerObjective::makeFromLikelihood(likelihood, _prior);
    _parameters[ndarray::view(0, nonlinearDim)] = initialNonlinear;
    _parameters[ndarray::view(nonlinearDim, parameterDim)] = initialAmplitudes;
    Optimizer optimizer(objective, _parameters, _ctrl.optimizer);
    int state = 0;
    if (_ctrl.doRecordHistory) {
        result._history = afw::table::BaseCatalog(_historySchema);
        state = optimizer.run(*_historyRecorder, result._history);
    } else {
        state = optimizer.run();
    }
    if (state & Optimizer::CONVERGED_TR_SMALL) {
        result._flags.set(Result::TR_SMALL);
    } else if ((state & Optimizer::FAILED_MAX_INNER_ITERATIONS)
               || (state & Optimizer::FAILED_MAX_OUTER_ITERATIONS)) {
        result._flags.set(Result::MAX_ITERATIONS);
    }
    result._nonlinear = optimizer.getParameters()[ndarray::view(0, nonlinearDim)];
    result._amplitudes = optimizer.getParameters()[ndarray::view(nonlinearDim, parameterDim)];
    result._hessian = optimizer.getHessian();
    result._objectiveValue = optimizer.getObjectiveValue();
    result._ellipses = _model->writeEllipses(result._nonlinear, fixed);
    for (
        Model::EllipseVector::iterator iter = result._ellipses.begin();
        iter != result._ellipses.end();
        ++iter
    ) {
        iter->transform(result._fitSysToMeasSys.geometric).inPlace();
    }
    result._flux = result._amplitudes.asEigen().sum() * result._fitSysToMeasSys.flux;
    result._fluxSigma = std::numeric_limits<Scalar>::quiet_NaN(); // tricky; will do later
    return result;
}

}}} // namespace lsst::meas::multifit
