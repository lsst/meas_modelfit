// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_SOURCE_MEASUREMENT_H
#define LSST_MEAS_MULTIFIT_SOURCE_MEASUREMENT_H
#include "boost/cstdint.hpp"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/afw/math/shapelets/constants.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include <Eigen/Core>

namespace lsst {
namespace meas {
namespace multifit {

/**
 *  @brief Options for configuring SourceMeasurement.
 *
 *  Should generally be referred to as SourceMeasurement::Options; it's not an inner class
 *  only because SWIG chokes on them.
 */
struct SourceMeasurementOptions {
    bool fitDeltaFunction;        ///< Whether to include a point-source component in the model.
    bool fitExponential;          ///< Whether to include an n=1 Sersic component in the model.
    bool fitDeVaucouleur;         ///< Whether to include an n=4 Sersic component in the model.
    int shapeletOrder;            ///< Order of additional shapelet expansion (none if < 0).
    int psfShapeletOrder;         ///< Shapelet order for PSF model approximation.
    int nGrowFp;                  ///< Number of pixels to grow footprints by.
    bool usePixelWeights;         ///< Whether to use per-pixel inverse variances as weights.
    double ellipticityStepSize;   ///< Range of ellipticities in brute-force grid.
    int ellipticityStepCount;     ///< Number of ellipticity points to test in each direction.
    double radiusMaxFactor;       ///< Sets upper limit of radius points in brute-force grid.
    double radiusMinFactor;       ///< Sets lower limit of radius points in brute-force grid.
    int radiusStepCount;          ///< Number of radius points to test.
    std::vector<std::string> maskPlaneNames; ///< Image mask bits to remove from footprint.
};

/**
 *  Guts of the ShapeletModelPhotometry measurement algorithm, with intermediate results
 *  stored as data members so we can inspect them and/or save them.
 *
 *  The model is optimized on a very sparse grid in e1, e2, and radius.  The grid points
 *  are set with the following procedure:
 *   - Define a pair of max and min radius ellipses by scaling the adaptive moments
 *     ellipse by radiusMaxFactor and radiusMinFactor.
 *   - Subtract the PSF moments from the max and min radius ellipse moments.  If either
 *     or both results in negative moments, set it to a zero-ellipticity, zero-radius
 *     ellipse.
 *   - Linearly interpolate radiusStepCount points between the minimum and maximum
 *     ellipses (note that we interpolate all three ellipse parameters, not just radius).
 *   - At each point, perturb the ellipticity in steps of ellipticityStepSize in each
 *     dimension and each direction (so we actually test (ellipticityStepCount * 2 + 1)^2
 *     points for each radius).  If the radius is zero, only zero ellipticity is tested.
 */
class SourceMeasurement {
public:
    typedef SourceMeasurementOptions Options;

    SourceMeasurement(Options const & options);

    static CompoundShapeletModelBasis::Ptr loadBasis(std::string const & name);

    static CompoundShapeletModelBasis::Ptr getExponentialBasis();
    static CompoundShapeletModelBasis::Ptr getDeVaucouleurBasis();

    void addObjectsToDefinition(
        Definition & definition, lsst::afw::geom::ellipses::Ellipse const & ellipse
    ) const;

    static Options readPolicy(lsst::pex::policy::Policy const & policy);

    static int computeCoefficientSize(Options const & options) {
        return options.fitDeltaFunction + options.fitExponential + options.fitDeVaucouleur
            + ((options.shapeletOrder < 0) ? 0 : afw::math::shapelets::computeSize(options.shapeletOrder));
    }

    static lsst::afw::geom::ellipses::Ellipse makeEllipse(
        lsst::afw::detection::Source const & source,
        lsst::afw::detection::Footprint const & fp
    );

    template <typename ExposureT>
    int measure(
        PTR(ExposureT) exp,
        PTR(lsst::afw::detection::Source) src
    ) {
        CONST_PTR(ExposureT) const_exp(exp);
        CONST_PTR(lsst::afw::detection::Source) const_src(src);
        return measure(const_exp, const_src);
    }

#ifndef SWIG
    template <typename ExposureT>
    int measure(
        CONST_PTR(ExposureT) exp,
        CONST_PTR(lsst::afw::detection::Source)
    );
#endif

    Evaluator::Ptr getEvaluator() const { return _evaluator; }
    lsst::afw::detection::Footprint::Ptr getFootprint() const { return _fp; }
    lsst::ndarray::Array<double const,1,1> getParameters() const { return _parameters; }
    lsst::ndarray::Array<double const,1,1> getCoefficients() const { return _coefficients; }
    lsst::ndarray::Array<double const,2,2> getCovariance() const { return _covariance; }
    lsst::ndarray::Array<double const,1,1> getIntegration() const { return _integration; }

    int getCoefficientSize() const { return _coefficients.getSize<0>(); }

    double getDeltaFunctionCoefficients() const {
        return _coefficients[getCoefficientOffset(DELTAFUNCTION_ID)];
    }
    double getExponentialCoefficients() const {
        return _coefficients[getCoefficientOffset(EXPONENTIAL_ID)];
    }
    double getDeVaucouleurCoefficients() const {
        return _coefficients[getCoefficientOffset(DEVAUCOULEUR_ID)];
    }

    lsst::ndarray::Array<double const,1,1> getShapeletCoefficients() const {
        int offset = getCoefficientOffset(SHAPELET_ID);
        return _coefficients[
            ndarray::view(offset, offset + afw::math::shapelets::computeSize(_options.shapeletOrder))
        ];
    }

    double getFlux() const { return _flux; }
    double getFluxErr() const { return _fluxErr; }
    Ellipse const & getEllipse() const { return _ellipse; }
    boost::int64_t getStatus() const { return _status; }

    CompoundShapeletModelBasis::Ptr getShapeletBasis() const { return _shapeletBasis; } 
    lsst::afw::image::MaskPixel getBitmask() const { return _bitmask; }
    Options const & getOptions() const { return _options; }

private:

    int getCoefficientOffset(ID id) const {
        return grid::find(_evaluator->getGrid()->objects, id).getCoefficientOffset();
    }

    void optimize(Ellipse const & initialEllipse);
    void solve(double e1, double e2, double radius, double & best);

    Options _options;
    lsst::afw::image::MaskPixel _bitmask;
    CompoundShapeletModelBasis::Ptr _shapeletBasis;

    double _flux, _fluxErr;
    Ellipse _ellipse;
    Evaluator::Ptr _evaluator;
    lsst::afw::detection::Footprint::Ptr _fp;
    ndarray::Array<double,1,1> _parameters;
    ndarray::Array<double,1,1> _coefficients;
    ndarray::Array<double,2,2> _covariance;
    ndarray::Array<double,1,1> _integration;
    boost::int64_t _status;

    static ID const DELTAFUNCTION_ID=0;
    static ID const EXPONENTIAL_ID=1;
    static ID const DEVAUCOULEUR_ID=2;
    static ID const SHAPELET_ID=3;

};

}}}

#endif
