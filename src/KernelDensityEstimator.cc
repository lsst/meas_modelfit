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

#include "ndarray.h"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/KernelDensityEstimator.h"
#include "lsst/meas/multifit/SampleSet.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

// maximum extent to evaluate KDE kernels, in units of sigma
double const KDE_MAX_RADIUS = 4.0;

// Gaussian FWHM in units of sigma
double const FWHM_FACTOR = 2.0 * std::sqrt(2.0 * std::log(2.0));

// Constant used in Gaussian normalization
double const SQRT_2PI = std::sqrt(2.0 * M_PI);

double computeStepAndCheck(double min, double max, int count, double fwhm) {
    if (min >= max) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Minimum range of histogram (%f) is above maximum (%f)") % min % max).str()
        );
    }
    if (count < 1) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Number of histogram bins (%d) must be greater than 0") % count).str()
        );
    }
    double step = (max - min) / count;
    if (!(fwhm > 2.0*step)) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("FWHM of kernel (%f) is not band limited (%f)")
             % fwhm % (2.0*step)).str()
        );
    }
    return step;
}

} // anonymous

KernelDensityEstimatorControl::KernelDensityEstimatorControl(
    int index, double min, double max, int count, double fwhm
) :
    _index(index), _count(count), _min(min), _max(max), _fwhm(fwhm),
    _step(computeStepAndCheck(min, max, count, fwhm))
{}

KernelDensityEstimatorControl::KernelDensityEstimatorControl(
    int index, double min, double max, int count
) :
    _index(index), _count(count), _min(min), _max(max), _fwhm(3.0 * (_max - _min) / _count ),
    _step(computeStepAndCheck(min, max, count, _fwhm))
{}

double KernelDensityEstimatorControl::getSigma() const {
    return _fwhm / FWHM_FACTOR;
}

void KernelDensityEstimatorControl::setMin(double min) {
    _step = computeStepAndCheck(min, _max, _count, _fwhm);
    _min = min;
}

void KernelDensityEstimatorControl::setMax(double max) {
    _step = computeStepAndCheck(_min, max, _count, _fwhm);
    _max = max;
}

void KernelDensityEstimatorControl::setCount(int count) {
    _step = computeStepAndCheck(_min, _max, count, _fwhm);
    _count = count;
}

void KernelDensityEstimatorControl::setSigma(double sigma) {
    setFwhm(sigma * FWHM_FACTOR);
}

void KernelDensityEstimatorControl::setFwhm(double fwhm) {
    _step = computeStepAndCheck(_min, _max, _count, fwhm);
    _fwhm = fwhm;
}

ndarray::Array<double,1,1> KernelDensityEstimatorControl::getCenters() const {
    ndarray::Array<double,1,1> result = ndarray::allocate(_count);
    double position = _min + 0.5 * _step;
    for (ndarray::Array<double,1,1>::Iterator i = result.begin(); i != result.end(); ++i, position += _step) {
        *i = position;
    }
    return result;
}

ndarray::Array<double,1,1> KernelDensityEstimatorControl::getEdges() const {
    ndarray::Array<double,1,1> result = ndarray::allocate(_count + 1);
    double position = _min;
    for (ndarray::Array<double,1,1>::Iterator i = result.begin(); i != result.end(); ++i, position += _step) {
        *i = position;
    }
    result[_count] = _max;  // better to make final step subject to round-off error than final edge
    return result;
}

void KernelDensityEstimatorControl::apply(
    ndarray::Array<double,1,1> const & output, double position, double weight
) const {
    double const sigma = getSigma();
    int const iCenter = int(std::floor((position - _min)/_step));
    int const iRadius = int(std::ceil(sigma * KDE_MAX_RADIUS));
    int const iMin = std::max(iCenter - iRadius, 0);
    int const iMax = std::min(iCenter + iRadius, _count - 1);
    weight /= (sigma * SQRT_2PI);
    for (int i = iMin; i <= iMax; ++i) {
        double x = (_min + _step * i - position) / sigma;
        output[i] += weight * std::exp(-0.5*x*x);
    }
}

}}} // namespace lsst::meas::multifit
