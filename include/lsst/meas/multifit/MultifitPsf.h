#ifndef LSST_MEAS_MULTIFIT_MULTIFITPSF_H
#define LSST_MEAS_MULTIFIT_MULTIFITPSF_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>
#include <complex>

#include "ndarray.hpp"

namespace lsst {
namespace meas {
namespace multifit {

typedef Eigen::Vector2d Coordinate;

class Kernel {
public:
    typedef boost::shared_ptr<Kernel> Ptr;
    typedef boost::shared_ptr<const Kernel> ConstPtr;

    // Should add to image rather than overwrite it.
    virtual void getImage(
        ndarray::ArrayRef<double,2,1> const & image, 
        Coordinate const & center
    ) const = 0;

    // Will assume image is zeros.
    virtual void getCentroidDerivative(
        ndarray::ArrayRef<double,2,1> const & dx,
        ndarray::ArrayRef<double,2,1> const & dy,
        Coordinate const & center
    ) const = 0;


    // Should add to output rather than overwrite it.
    virtual void convolve(
        ndarray::ArrayRef<double,2,2> const & input,
        ndarray::ArrayRef<double,2,2> const & output,
        int oversampling = 1
    ) const;
    
    // Should add to output rather than overwrite it.
    virtual void fourierConvolve(
        ndarray::ArrayRef<std::complex<double>,2,2> const & input,
        ndarray::ArrayRef<double,2,2> const & output,
        int oversampling = 1
    ) const;

    virtual ~Kernel() {}

protected:
    virtual void getArray(ndarray::ArrayRef<double,2,2> const & image, int oversampling) const = 0;
    virtual void getFourierArray(ndarray::ArrayRef<std::complex<double>,2,2> const & image,
                                 int oversampling) const;
    void operator=(Kernel const & other) {}
};

class Psf {
public:

    typedef boost::shared_ptr<Psf> Ptr;
    typedef boost::shared_ptr<const Psf> ConstPtr;
    
    virtual int getBasisSize() const = 0;
    
    virtual Kernel::Ptr getKernel() const;
    
    virtual Kernel::Ptr getBasisKernel(int order) const = 0;
    
    virtual Eigen::VectorXd getCoefficients() const = 0;
    
    virtual Eigen::MatrixXd getCovariance() const = 0;

    virtual ~Psf() {}

protected:
    void operator=(Psf const & other) {}
};

}

#endif // !MODELING_Psf_hpp_INCLUDED
