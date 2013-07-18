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

#ifndef LSST_MEAS_MULTIFIT_KernelDensityEstimator_h_INCLUDED
#define LSST_MEAS_MULTIFIT_KernelDensityEstimator_h_INCLUDED

#include "ndarray_fwd.h"

namespace lsst { namespace meas { namespace multifit {

class SampleSet;

/**
 *  @brief Control object used along with SampleSet::computeDensity to create a kernel-density estimate
 *         (roughly, a smoothed histogram) of a dimension of the samples.
 *
 *  See e.g. https://en.wikipedia.org/wiki/Kernel_density_estimation
 */
class KernelDensityEstimatorControl {
public:

    /// Basic constructor; see corresponding accessors for argument descriptions.
    KernelDensityEstimatorControl(int index, double min, double max, int count, double fwhm);

    /// Construct while automatically setting the FWHM to 3 times the step size.
    KernelDensityEstimatorControl(int index, double min, double max, int count);

    /// Return the index of the parameter to which this histogram range applies
    int getIndex() const { return _index; }

    /// Return the lower edge of the first bin the histogram
    double getMin() const { return _min; }

    /// Return the upper edge of the last bin the histogram
    double getMax() const { return _max; }

    /// Return the distance between grid points
    double getStep() const { return _step; }

    /// Return the number of grid points
    int getCount() const { return _count; }

    /// Return the width of the Gaussian kernel
    double getSigma() const;

    /// Return the width of the Gaussian kernel
    double getFwhm() const { return _fwhm; }

    /// Set the index of the parameter to which this histogram range applies
    void setIndex(int index) { _index = index; }

    /// Set the lower edge of the first bin, modifying bin size to keep bin count and max fixed
    void setMin(double min);

    /// Set the upper edge of the last bin, modifying bin size to keep bin count and min fixed
    void setMax(double max);

    /// Set the number of grid points, modifying bin size to keep min and max fixed
    void setCount(int count);

    /// Set the width of the Gaussian kernel; must always be wide enough for a Nyquist-sampled grid
    void setSigma(double sigma);

    /// Set the width of the Gaussian kernel; must always be wide enough for a Nyquist-sampled grid
    void setFwhm(double fwhm);

    /// Return the centers of the grid pixels (size == getCount())
    ndarray::Array<double,1,1> getCenters() const;

    /// Return the edges of the grid pixels (size == getCount() + 1)
    ndarray::Array<double,1,1> getEdges() const;

protected:

    friend class SampleSet;

    void apply(ndarray::Array<double,1,1> const & output, double position, double weight) const;

    int _index;
    int _count;
    double _min;
    double _max;
    double _fwhm;
    double _step;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_DensityEstimators_h_INCLUDED
