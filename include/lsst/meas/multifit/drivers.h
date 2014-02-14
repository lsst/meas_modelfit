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

struct OptimizerFitControl {

    LSST_NESTED_CONTROL_FIELD(
        likelihood, lsst.meas.multifit, ProjectedLikelihoodControl,
        "configuration for how the model is evaluated and residuals are weighted"
    );

    LSST_NESTED_CONTROL_FIELD(
        optimizer, lsst.meas.multifit, OptimizerControl,
        "configuration for how the objective surface is explored"
    );

};

class OptimizerFit {
public:

    typedef OptimizerFitControl Control;

    OptimizerFit(
        PTR(Model) model, PTR(Prior) prior,
        PTR(afw::coord::Coord) position,
        UnitSystem const & fitSys, UnitSystem const & measSys,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar const,1,1> const & fixed,
        Control const & ctrl
    );

    OptimizerFit(OptimizerFit const & other);

    OptimizerFit & operator=(OptimizerFit const & other);

    void swap(OptimizerFit & other);

    void run(
        afw::image::Exposure<Pixel> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        bool doRecordHistory=true
    );

    //@{
    /**
     *  @brief Accessors for the Model to be fit to the data.
     *
     *  While the Model can be changed after the fitter has been constructed, the new Model must have
     *  the same amplitude dimension as the last one, and if it has different nonlinear or fixed
     *  dimensions, it must have the same number of ellipses.  These constraints allow the parameters
     *  to be carried over from one model to the other, so one can e.g. start by fitting a simple,
     *  approximate model with some parameters fixed, before fine-tuning with a more complex model.
     *
     *  If setModel() is called with doForceEllipseConversion=true, the parameters will be converted
     *  through conversion to ellipses even if the number of nonlinear and fixed parameters is unchanged.
     */
    PTR(Model) getModel() const { return _model; }
    void setModel(PTR(Model) model, bool doForceEllipseConversion=false);
    //@}

    //@{
    /// Accessors for the Bayesian prior used to augment the likelihood and regularize the fit.
    PTR(Prior) getPrior() const { return _prior; }
    void setPrior(PTR(Prior) prior);
    //@}

    //@{
    // Accessors for the configuration of the fitter
    Control const & getControl() const { return _ctrl; }
    Control & getControl() { return _ctrl; }
    void setControl(Control const & ctrl) { _ctrl = ctrl; }
    //@}

    //@{
    /**
     *  @brief Accessors for the position of the object being fit.
     *
     *  This position is used when computing the FitSysToMeasSys transform (we need a position
     *  in order to evaluate the local linear approximations to the more general transforms),
     *  and setting it will cause that transform to be updated.
     */
    PTR(afw::coord::Coord) getPosition() const { return _position; }
    void setPosition(PTR(afw::coord::Coord) position);
    //@}

    /// Return value of the objective function (likelihood or unnormalized posterior) at the best-fit point.
    Scalar getObjectiveValue() const { return _objectiveValue; }

    /// Return the bitflags that describe the state of the optimizer (see Optimizer::StateFlags).
    int getOptimizerState() const { return _optimizerState; }

    /// Return the current nonlinear parameters, in the "FitSys" UnitSystem.
    ndarray::Array<Scalar const,1,1> getNonlinear() const { return _nonlinear; }

    /// Return the current amplitude parameters, in the "FitSys" UnitSystem.
    ndarray::Array<Scalar const,1,1> getAmplitudes() const { return _amplitudes; }

    /// Return the fixed parameters (values specific to this object, but held fixed during all fits).
    ndarray::Array<Scalar const,1,1> getFixed() const { return _fixed; }

    /// Return the flux in the "MeasSys" UnitSystem.
    Scalar getFlux() const;

    /// Return the 1-sigma uncertainty on the flux in the "MeasSys" UnitSystem.
    Scalar getFluxSigma() const;

    /// Return the nth ellipse of the model, in the "MeasSys" UnitSystem.
    afw::geom::ellipses::Ellipse getEllipse(int n=0) const;

    /// Return the UnitSystem that is used to define the parameters in the optimizer
    UnitSystem getFitSys() const { return _fitSys; }

    /// Return the UnitSystem in which measured quantities (i.e. flux, ellipses) should be reported
    UnitSystem getMeasSys() const { return _measSys; }

    /// Return the transforms from FitSys to MeasSys at the position of the object being fit.
    LocalUnitTransform getFitSysToMeasSys() const { return _fitSysToMeasSys; }

    /// Return a catalog that traces the optimizer path (if doRecordHistory=true on the last call to run())
    afw::table::BaseCatalog getHistory() const { return _history; }

private:

    void _ensureMeasQuantities() const;

    PTR(Model) _model;
    PTR(Prior) _prior;
    PTR(afw::coord::Coord) _position;
    mutable bool _hasMeasQuantities;
    int _optimizerState;
    Scalar _objectiveValue;
    mutable Scalar _flux;
    mutable Scalar _fluxSigma;
    UnitSystem _fitSys;
    UnitSystem _measSys;
    LocalUnitTransform _fitSysToMeasSys;
    mutable Model::EllipseVector _ellipses;
    ndarray::Array<Scalar,1,1> _parameters;
    ndarray::Array<Scalar,1,1> _nonlinear;
    ndarray::Array<Scalar,1,1> _amplitudes;
    ndarray::Array<Scalar,1,1> _fixed;
    ndarray::Array<Scalar,2,2> _hessian;
    PTR(OptimizerHistoryRecorder) _historyRecorder;
    afw::table::BaseCatalog _history;
    Control _ctrl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_drivers_h_INCLUDED
