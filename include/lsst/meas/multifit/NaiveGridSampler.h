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

#ifndef LSST_MEAS_MULTIFIT_NaiveGridSampler_h_INCLUDED
#define LSST_MEAS_MULTIFIT_NaiveGridSampler_h_INCLUDED

#include "lsst/meas/multifit/BaseSampler.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief SamplerState class that evaluates on a simple grid in (radius, e1, e2).
 *
 *  For each radius from zero to maxRadius, we evaluate each grid point (e_1, e_2) with spacing
 *  ellipticityStepSize, starting from (0, 0), such that |e| < maxEllipticity.
 */
class NaiveGridSampler : public BaseSampler {
public:

    typedef afw::geom::ellipses::Separable<
        afw::geom::ellipses::ReducedShear,
        afw::geom::ellipses::TraceRadius
    > EllipseCore;

    /**
     *  @brief Compute the total number of samples on the grid.
     */
    static int computeSampleSetSize(
        int nRadiusSteps,
        double ellipticityStepSize,
        double maxEllipticity
    );

    /**
     *  @brief Construct with grid size parameters.
     *
     *  @param[in] center               Center point of the elliptical model for all samples
     *  @param[in] nRadiusSteps         Number of steps in radius on the grid
     *  @param[in] ellipticityStepSize  Size of step in e1 and e2
     *  @param[in] maxRadius            Maximum radius value
     *  @param[in] maxEllipticity       Maximum ellipticity magnitude value
     *
     *  We use the number of radius steps because we expect the maximum radius to be different
     *  for each object, but we may want to make the total number of grid points constant.  We
     *  do not expect the maximum ellipticity to be different for each object, so we use a
     *  more natural step size parameter to set the grid there.
     */
    NaiveGridSampler(
        afw::geom::Point2D const & center,
        int nRadiusSteps,
        double ellipticityStepSize,
        double maxRadius,
        double maxEllipticity
    );

    /**
     *  @brief Generate and evaluate samples on the grid
     *
     *  See the NaiveGridSampler class documentation for the details of the grid.
     */
    virtual SampleSet run(Objective const & objective) const;

private:
    afw::geom::Point2D _center;
    int _nRadiusSteps;
    double _ellipticityStepSize;
    double _maxRadius;
    double _maxEllipticity;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_NaiveGridSampler_h_INCLUDED
