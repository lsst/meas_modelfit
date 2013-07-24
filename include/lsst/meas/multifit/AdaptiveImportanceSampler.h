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

#ifndef LSST_MEAS_MULTIFIT_AdaptiveImportanceSampler_h_INCLUDED
#define LSST_MEAS_MULTIFIT_AdaptiveImportanceSampler_h_INCLUDED

#include "lsst/meas/multifit/BaseSampler.h"

namespace lsst { namespace meas { namespace multifit {

class ImportanceSamplerControl {
public:
    LSST_CONTROL_FIELD(nSamples, int, "Number of Monte Carlo samples to draw");
    LSST_CONTROL_FIELD(nUpdateSteps, int, "Number of Expectation-Maximization update iterations");
    LSST_CONTROL_FIELD(
        keepPrevious, bool,
        "Whether to keep previous samples and append to them, adjusting their probability accordingly"
    );
};

class AdaptiveImportanceSamplerControl {
public:
    LSST_CONTROL_FIELD(
        degreesOfFreedom, double,
        "Number of degrees of freedom in Student's T mixture components (inf=Gaussian)"
    );

    LSST_NESTED_CONTROL_FIELD(
        iteration1, lsst.meas.multifit, ImportanceSamplerControl,
        "Config parameters for first importance sampling iteration"
    );
    LSST_NESTED_CONTROL_FIELD(
        iteration2, lsst.meas.multifit, ImportanceSamplerControl,
        "Config parameters for second importance sampling iteration"
    );
    LSST_NESTED_CONTROL_FIELD(
        iteration3, lsst.meas.multifit, ImportanceSamplerControl,
        "Config parameters for third importance sampling iteration"
    );

};


/**
 *  @brief SamplerState class that evaluates on a simple grid in (radius, e1, e2).
 *
 *  For each radius from zero to maxRadius, we evaluate each grid point (e_1, e_2) with spacing
 *  ellipticityStepSize, starting from (0, 0), such that |e| < maxEllipticity.
 */
class AdaptiveImportanceSampler : public BaseSampler {
public:

    /**
     *  @brief Generate and evaluate samples on the grid
     *
     *  See the NaiveGridSampler class documentation for the details of the grid.
     */
    virtual SampleSet run(Objective const & objective) const;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_AdaptiveImportanceSampler_h_INCLUDED
