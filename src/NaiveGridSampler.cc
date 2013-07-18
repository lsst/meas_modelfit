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
#include "lsst/meas/multifit/NaiveGridSampler.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

// Function to iterate over all ellipticities in the grid, to we ensure that we use
// exactly the same grid everywhere
template <typename Functor>
void forEachEllipticity(Functor & f, double maxEllipticity, double ellipticityStepSize) {
    // we compare squared values so we can make e1 and e2 behave more similarly in the precence of
    // round-off error
    double const e1maxSqr = maxEllipticity * maxEllipticity;
    {
        f(0.0, 0.0);
        double e2maxSqr = e1maxSqr;
        for (double e2 = ellipticityStepSize; e2*e2 < e2maxSqr; e2 += ellipticityStepSize) {
            f(0.0, -e2);
            f(0.0, e2);
        }
    }
    for (double e1 = ellipticityStepSize; e1*e1 < e1maxSqr; e1 += ellipticityStepSize) {
        f(-e1, 0.0);
        f(e1, 0.0);
        double e2maxSqr = e1maxSqr - e1 * e1;
        for (double e2 = ellipticityStepSize; e2*e2 < e2maxSqr; e2 += ellipticityStepSize) {
            f(-e1, -e2);
            f(e1, -e2);
            f(-e1, e2);
            f(e1, e2);
        }
    }
}

struct CountSamples {
    void operator()(double e1, double e2) { ++n; }

    CountSamples() : n(0) {}

    int n;
};

struct AddConstantSample {

    void operator()(double e1, double e2) {
        (*parameters)[NaiveGridSampler::EllipseCore::E1] = e1;
        (*parameters)[NaiveGridSampler::EllipseCore::E2] = e2;
        samples->add(*joint, proposal, *parameters);
    }

    AddConstantSample(
        LogGaussian * joint_, samples::Scalar proposal_, samples::Vector * parameters_, SampleSet * samples_
    ) : joint(joint_), proposal(proposal_), parameters(parameters_), samples(samples_) {}

    LogGaussian * joint;
    samples::Scalar proposal;
    samples::Vector * parameters;
    SampleSet * samples;
};

struct AddNewSample {

    void operator()(double e1, double e2) {
        NaiveGridSampler::EllipseCore & ec
            = static_cast<NaiveGridSampler::EllipseCore &>(ellipse->getCore());
        ec.setE1(e1);
        ec.setE2(e2);
        (*parameters)[NaiveGridSampler::EllipseCore::E1] = e1;
        (*parameters)[NaiveGridSampler::EllipseCore::E2] = e2;
        *joint = objective->evaluate(*ellipse);
        samples->add(*joint, proposal, *parameters);
    }

    AddNewSample(
        LogGaussian * joint_, samples::Scalar proposal_, samples::Vector * parameters_,
        SampleSet * samples_,
        afw::geom::ellipses::Ellipse * ellipse_,
        Objective const * objective_
    ) : joint(joint_), proposal(proposal_), parameters(parameters_),
        samples(samples_), ellipse(ellipse_), objective(objective_) {}

    LogGaussian * joint;
    samples::Scalar proposal;
    samples::Vector * parameters;
    SampleSet * samples;
    afw::geom::ellipses::Ellipse * ellipse;
    Objective const * objective;
};

} // anonymous


int NaiveGridSampler::computeSampleSetSize(
    int nRadiusSteps,
    double ellipticityStepSize,
    double maxEllipticity
) {
    CountSamples f;
    forEachEllipticity(f, maxEllipticity, ellipticityStepSize);
    return f.n * nRadiusSteps;
}

NaiveGridSampler::NaiveGridSampler(
    afw::geom::Point2D const & center,
    int nRadiusSteps,
    double ellipticityStepSize,
    double maxRadius,
    double maxEllipticity
) : _center(center),
    _nRadiusSteps(nRadiusSteps),
    _ellipticityStepSize(ellipticityStepSize),
    _maxRadius(maxRadius),
    _maxEllipticity(maxEllipticity)
{
    if (_nRadiusSteps < 2) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "nRadiusSteps must be >= 2"
        );
    }
    if (_maxEllipticity < 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "maxEllipticity must be >= 0.0"
        );
    }
    if (_ellipticityStepSize < 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "ellipticityStepSize must be >= 0.0"
        );
    }
    if (!(_maxRadius >= 0.0)) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            "maxRadius must be >= 0.0"
        );
    }
}

SampleSet NaiveGridSampler::run(Objective const & objective) const {
    int nonlinearDim = 3;
    int linearDim = objective.getLinearDim();
    afw::geom::ellipses::Ellipse ellipse = afw::geom::ellipses::Ellipse(EllipseCore(), _center);
    SampleSet samples(nonlinearDim, linearDim, ellipse.getCore().getName());
    samples.setDataSquaredNorm(objective.getDataSquaredNorm());
    double proposal = std::log(_maxRadius * _maxEllipticity * _maxEllipticity * M_PI);
    ellipse.getCore().scale(0.0);
    samples::Vector parameters = ellipse.getCore().getParameterVector().head(nonlinearDim);
    LogGaussian joint = objective.evaluate(ellipse);
    // At r=0, we only evaluate e1=0, e2=0 and use that likelihood for all ellipticities.
    // It's a waste of space, but it makes analysis easier, and at least it's not a waste
    // of cycles.
    {
        AddConstantSample f(&joint, proposal, &parameters, &samples);
        // the loop over ellipticities happen inside this function, which calls f repeatedly.
        forEachEllipticity(f, _maxEllipticity, _ellipticityStepSize);
    }
    double radiusStepSize = _maxRadius / (_nRadiusSteps - 1);
    for (int i = 1; i < _nRadiusSteps; ++i) {
        parameters[EllipseCore::RADIUS] = radiusStepSize * i;
        static_cast<EllipseCore &>(ellipse.getCore()).setRadius(parameters[EllipseCore::RADIUS]);
        AddNewSample f(&joint, proposal, &parameters, &samples, &ellipse, &objective);
        forEachEllipticity(f, _maxEllipticity, _ellipticityStepSize);
    }
    return samples;
}

}}} // namespace lsst::meas::multifit
