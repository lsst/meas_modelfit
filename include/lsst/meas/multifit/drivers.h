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

#ifndef LSST_MEAS_MULTIFIT_drivers_h_INCLUDED
#define LSST_MEAS_MULTIFIT_drivers_h_INCLUDED

#include <bitset>

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/ProjectedLikelihood.h"
#include "lsst/meas/multifit/optimizer.h"

namespace lsst { namespace meas { namespace multifit {

struct OptimizerFitterControl {

    LSST_NESTED_CONTROL_FIELD(
        likelihood, lsst.meas.multifit, ProjectedLikelihoodControl,
        "configuration for how the model is evaluated and residuals are weighted"
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.multifit, OptimizerControl,
        "configuration for how the objective surface is explored"
    );

    LSST_CONTROL_FIELD(
        doRecordHistory, bool,
        "whether to record the steps the optimizer takes for debugging purposes"
    );

};

class OptimizerFitterResult {
public:

    enum FlagBit {
        MAX_ITERATIONS=0,
        TR_SMALL,
        N_FLAGS
    };

    OptimizerFitterResult(
        PTR(Model) model, afw::coord::Coord const & position,
        UnitSystem const & fitSys, UnitSystem const & measSys
    ) :
        _model(model), _fitSys(fitSys), _measSys(measSys), _fitSysToMeasSys(position, _fitSys, _measSys)
    {}

    Scalar getObjectiveValue() const { return _objectiveValue; }

    bool getFlag(FlagBit bit) const { return _flags[bit]; }

    ndarray::Array<Scalar const,1,1> getNonlinear() const { return _nonlinear; }

    ndarray::Array<Scalar const,1,1> getAmplitudes() const { return _amplitudes; }

    ndarray::Array<Scalar const,1,1> getFixed() const { return _fixed; }

    UnitSystem getFitSys() const { return _fitSys; }

    afw::table::BaseCatalog getHistory() const { return _history; }

    Scalar getFlux() const { return _flux; }

    Scalar getFluxSigma() const { return _fluxSigma; }

    afw::geom::ellipses::Ellipse getEllipse(int i=0) const { return _ellipses[i]; }

protected:

    friend class OptimizerFitter;

    PTR(Model) _model;
    Scalar _objectiveValue;
    Scalar _flux;
    Scalar _fluxSigma;
    UnitSystem _fitSys;
    UnitSystem _measSys;
    LocalUnitTransform _fitSysToMeasSys;
    Model::EllipseVector _ellipses;
    ndarray::Array<Scalar const,1,1> _nonlinear;
    ndarray::Array<Scalar const,1,1> _amplitudes;
    ndarray::Array<Scalar const,1,1> _fixed;
    ndarray::Array<Scalar const,2,2> _hessian;
    std::bitset<N_FLAGS> _flags;
    afw::table::BaseCatalog _history;
};

class OptimizerFitter {
public:

    typedef OptimizerFitterControl Control;
    typedef OptimizerFitterResult Result;

    OptimizerFitter(PTR(Model) model, PTR(Prior) prior, Control const & ctrl);

    ndarray::Array<Scalar const,1,1> getDefaultInitialAmplitudes() const;

    Result apply(
        afw::image::Exposure<Pixel> const & data,
        afw::detection::Footprint const & footprint,
        afw::coord::Coord const & position,
        UnitSystem const & fitSys,
        UnitSystem const & measSys,
        ndarray::Array<Scalar const,1,1> const & initialNonlinear,
        ndarray::Array<Scalar const,1,1> const & initialAmplitudes,
        ndarray::Array<Scalar const,1,1> const & fixed,
        shapelet::MultiShapeletFunction const & psf
    ) const;

private:
    PTR(Model) _model;
    PTR(Prior) _prior;
    afw::table::Schema _historySchema;
    PTR(OptimizerHistoryRecorder) _historyRecorder;
    ndarray::Array<Scalar,1,1> _parameters;
    Control _ctrl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_drivers_h_INCLUDED
