#ifndef LSST_MEAS_MULTIFIT_CALIBRATED_EXPOSURE
#define LSST_MEAS_MULTIFIT_CALIBRATED_EXPOSURE

#include "boost/shared_ptr.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"

#include "ndarray/ndarray.hpp"

namespace lsst {
namespace meas {
namespace multifit {

class Psf {
public:

    virtual int getNumParam() const = 0;

    virtual int getOversampling() const = 0;

    virtual ndarray::Index<2> getImageShape() const = 0;
    virtual ndarray::Index<2> getFourierImageShape() const = 0;
    
    // shape is (param,y,x)
    virtual ndarray::ArrayRef<double,3,3> getBasis() const = 0;

    // shape is (param,ky,kx)
    virtual ndarray::ArrayRef<std::complex<double>,3,3> getFourierBasis() const = 0;
    
    // shape is (y,x); equivalent to matrix product of parameters with basis
    virtual void fillImageAt(double y, double x,
            ndarray::ArrayRef<double,2,2> const & output) const = 0;

    // shape is (ky,kx); equivalent to matrix product of parameters with Fourier basis
    virtual void fillFourierImageAt(Eigen::Vector2d const & local_pos,
                                    ndarray::ArrayRef<std::complex<double>,2,2> const & output) const = 0;

    // shape is (param)
    virtual void fillParametersAt(Eigen::Vector2d const & local_pos,
                                    ndarray::ArrayRef<double,1,1> const & params) const = 0;

    // shape is (param,param); Error matrix == inverse of covariance matrix
    virtual void fillErrorMatrixAt(Eigen::Vector2d const & local_pos,
                                   ndarray::ArrayRef<double,2,2> const & output) const = 0;

    virtual ~Psf() {}
};

class Background {
public:

    virtual ndarray::Index<2> getImageShape() const = 0;    

    virtual int getNumParam() const = 0;

    // shape is (y,x)
    virtual ndarray::ArrayRef<double,2,2> getImage() const = 0;

    // shape is (param,y,x)
    virtual ndarray::ArrayRef<double,3,3> getMatrix() const = 0;

    // shape is (param,param); Error matrix == inverse of covariance matrix
    virtual ndarray::ArrayRef<double,2,2> getErrorMatrix() const = 0;
    
    virtual ~Background() {}
};

class CalibratedExposure {
public:
    typedef boost::shared_ptr<CalibratedExposure> Ptr;

    virtual Eigen::Transform2d getTransform(Eigen::Vector2d const & global_pos) const = 0;
    virtual double getPhotometricFactorAt(Eigen::Vector2d const & local_pos) const = 0;

    virtual Psf * const & getPsf() const = 0;
    virtual int getNumPsfParam() const {
        Psf * psf = getPsf();
        if(psf == NULL)
            return 0;
        else return psf->getNumParam();
    }
    virtual Background * const getBackground() const = 0;
    virtual int getNumBackgroundParam() const {
        Background * background = getBackground();
        if(background == NULL)
            return 0;
        else return getBackground()->getNumParam();
    }

    virtual int getWidth() const = 0;
    virtual int getHeight() const = 0;
    virtual ndarray::Index<2> getImageShape() const = 0;

    virtual ndarray::ArrayRef<double,2,2> getImage() const = 0;
    virtual ndarray::ArrayRef<double,2,2> getWeight() const = 0;

    virtual ~CalibratedExposure() {}
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_CALIBRATED_EXPOSURE
