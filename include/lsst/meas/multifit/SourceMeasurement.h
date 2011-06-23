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
#include "lsst/meas/multifit/GaussNewtonOptimizer.h"
#include "lsst/meas/multifit/BruteForceSourceOptimizer.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include <Eigen/Core>

namespace lsst {
namespace meas {
namespace multifit {


class SourceMeasurement {
public:
    enum {
        NO_EXPOSURE=0x001, 
        NO_PSF=0x002, 
        NO_SOURCE=0x004, 
        NO_BASIS=0x008,
        NO_FOOTPRINT=0x010, 
        BAD_INITIAL_MOMENTS=0x020, 
        OPTIMIZER_FAILED=0x040,
        GALAXY_MODEL_FAILED=0x080,
        UNSAFE_INVERSION=0x100
    };
    SourceMeasurement(int basis, int psfShapeletOrder,
                      int nTestPoints, int nGrowFp, 
                      bool usePixelWeights, bool fitDeltaFunction,
                      bool isEllipticityActive, 
                      bool isRadiusActive, 
                      bool isPositionActive,
                      std::vector<std::string> const & maskPlaneNames);

    SourceMeasurement(ModelBasis::Ptr basis, int psfShapeletOrder,
                      int nTestPoints, int nGrowFp, 
                      bool usePixelWeights, bool fitDeltaFunction,
                      bool isEllipticityActive, 
                      bool isRadiusActive, 
                      bool isPositionActive,
                      lsst::afw::image::MaskPixel bitmask);


    static CompoundShapeletModelBasis::Ptr loadBasis(int basisSize);

    static CompoundShapeletModelBasis::Ptr loadBasis(std::string const & name);

    static lsst::afw::geom::ellipses::Ellipse makeEllipse(
        lsst::afw::detection::Source const & source,
        lsst::afw::detection::Footprint const & fp
    );
     
    template <typename ExposureT>
    int measure(
        CONST_PTR(ExposureT) exp,
        CONST_PTR(lsst::afw::detection::Source)
    ); 
   
    CONST_PTR(Evaluator) getEvaluator() const {return _evaluator;}
    CONST_PTR(lsst::afw::detection::Footprint) getFootprint() const {return _fp;}
    ndarray::Array<double const, 1 ,1> getParam() const{return _param;}
    ndarray::Array<double const, 1, 1> getCoeff() const{return _coeff;}
    ndarray::Array<double const, 1, 1> getBasisCoeff() const {return _coeff[ndarray::view(0,getBasisSize())];}
    double getFlux() const {return _flux;}
    double getFluxErr() const {return _fluxErr;}
    Ellipse::ConstPtr getEllipse() const {return _ellipse;}
    boost::int64_t getStatus() const {return _status;}


    //examine configuration bits
    CONST_PTR(ModelBasis) getBasis() const {return _basis;} 
    lsst::afw::image::MaskPixel getBitmask() const {return _bitmask;}
    int getNGrowFp() const {return _nGrowFp;}
    int getBasisSize() const {return _basis->getSize();}
    int getCoefficientSize() const {return _coeff.getSize<0>();}
    bool usePixelWeights() const {return _usePixelWeights;}
    bool fitDeltaFunction() const {return _fitDeltaFunction;}
    bool isEllipticityActive() const {return _isEllipticityActive;} 
    bool isRadiusActive() const {return _isRadiusActive;} 
    bool isPositionActive() const {return _isPositionActive;} 

private:
    bool _usePixelWeights;
    bool _fitDeltaFunction;
    bool _isEllipticityActive, _isRadiusActive, _isPositionActive;
    lsst::afw::image::MaskPixel _bitmask;
    int _nTestPoints;
    int _nGrowFp;
    int _psfShapeletOrder;
    ModelBasis::Ptr _basis;

    double _flux, _fluxErr;
    Ellipse::Ptr _ellipse;
    Evaluator::Ptr _evaluator;
    CONST_PTR(lsst::afw::detection::Footprint) _fp;
    ndarray::Array<double, 1, 1> _param, _coeff;
    boost::int64_t _status;

    static ID const GALAXY_ID=0;
    static ID const STAR_ID=1;
};

}}}

#endif
