// -*- lsst-c++ -*-
/**
 * @file
 * Collection of useful typedefs used throughout the lsst::meas::multifit
 * namespace
 */
#ifndef LSST_MEAS_MULTIFIT_CORE_H
#define LSST_MEAS_MULTIFIT_CORE_H

#include <Eigen/Core>

#include "lsst/afw/geom.h"
#include "lsst/afw/geom/deprecated.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace lsst {
namespace meas {
namespace multifit {

typedef double Pixel;
typedef double Parameter;

typedef Parameter * ParameterIterator;
typedef Parameter const * ParameterConstIterator;
typedef std::map<int, Parameter> ParameterMap;

typedef Eigen::aligned_allocator<char> Allocator;
typedef Eigen::VectorXd ParameterVector;
typedef Eigen::Map<Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> > MatrixMap;
typedef Eigen::Block<MatrixMap> MatrixMapBlock;
typedef Eigen::Map<Eigen::Matrix<Pixel, Eigen::Dynamic, 1> > VectorMap;

typedef lsst::meas::algorithms::PSF Psf;
typedef lsst::afw::math::Kernel Kernel;
typedef lsst::afw::image::Wcs Wcs;
typedef lsst::afw::detection::Footprint Footprint;

typedef boost::shared_ptr<Psf const> PsfConstPtr;
typedef boost::shared_ptr<Kernel const> KernelConstPtr;
typedef boost::shared_ptr<Wcs const> WcsConstPtr;
typedef boost::shared_ptr<Footprint const> FootprintConstPtr;

class ParameterRangeException: public lsst::pex::exceptions::RangeErrorException {
public:
    ParameterRangeException(
        char const * file, int line, char const * func,
        std::string const & message,
        ParameterMap outOfRangeParameterMap
    ) : lsst::pex::exceptions::RangeErrorException(file, line, func, message),
        _outOfRangeParameterMap(outOfRangeParameterMap) 
    {}
    
    virtual char const * getType(void) const throw() {
        return "lsst::meas::multifit::ParameterRangeException *";
    }
    virtual lsst::pex::exceptions::Exception* clone(void) const {
        return new ParameterRangeException(*this);
    }
    ParameterMap const & getOutOfRangeParameterMap() const {
        return _outOfRangeParameterMap;
    }
private:
    ParameterMap _outOfRangeParameterMap;
};

}}} // namespace lsst::meas::multifit

#endif
