#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include "boost/shared_ptr.hpp"

#include "Eigen/Core"
#include "ndarray/ndarray.hpp"
#include "lsst/meas/multifit/AffineTransform.h"

namespace lsst {
namespace meas {
namespace multifit {

struct ModelShape {
    int imageWidth, imageHeight;
    int linearSize;
    int nonlinearSize;
    int psfBasisSize;
    int transformSize;
}

typedef Eigen::Vector2d Coordinate;
typedef ndarray::ArrayRef<double, 1, 1> ParameterVector;
//dimensions are pixel height, pixel width
typedef ndarray::ArrayRef<double, 2, 2> ImageVector;
// dimensions are parameters, pixel height, pixel width
typedef ndarray::ArrayRef<double, 3, 3> DerivativeMatrix;

inline ImageVector getImageView(Eigen::VectorXd vector, 
        int const & height, int const & width
) {
    return ImageVector(vector.data(), ndarray::make_index(height, width));
}
inline ImageVector getImageView(Eigen::MatrixXd matrix, int row,
        int const & height, int const & width
) {
    return ImageVector(matrix.row(row).data(),
            ndarray::make_index(height, width));
}
inline DerivativeMatrix(Eigen::MatrixXd matrix, 
        int const & nParams, int const & height, int const width
) {
    return DerivativeMatrix(matrix.data(), 
            ndarray::make_index(nParameters, height, width));
} 

inline Eigen::Map<Eigen::VectorXd> extractEigenView(
        ParameterVector const & array
) {
    return array.core().vector();
}

inline Eigen::Map<Eigen::VectorXd> extractEigenView(
        ImageVector const & array
) {
    return array.core().vector();
}

inline Eigen::Map<Eigen::MatrixXd> extractEigenView(
        DerivativeMatrix const & array
) {
    return Eigen::Map<Eigen::MatrixXd>(
            array.data(),array.shape()[2]*array.shape()[1],array.shape()[0]
    );
}

class Model {
public:
    typedef boost::shared_ptr<Model> Ptr;
    typedef boost::shared_ptr<const Model> ConstPtr;

    virtual void setLinearParameters(Eigen::VectorXd const & parameters) = 0;
    virtual void setNonlinearParameters(Eigen::VectorXd const & parameters) = 0;
    virtual Eigen::VectorXd getLinearParameters() const = 0;  
    virtual Eigen::VectorXd getNonlinearParameters() const = 0;
    
    /**
     * Apply transform to the model "after" any existing transform.
     */
    virtual void addTransform(AffineTransform const & transform) = 0;
    virtual void setTransform(AffineTransform const & transform) = 0;
    virtual AffineTransform getTransform() const = 0;
    
    virtual Model * clone() const = 0;

    /**
     * Creates a convolved model
     * If the model already has a Psf, the old Psf is ignored.
     */
    virtual Model * convolve(Psf::ConstPtr psf) const = 0;

    Eigen::VectorXd const & computeParametrizedImage();
    Eigen::VectorXd const & computeConstantImage();
    Eigen::MatrixXd const & computeLinearMatrix();   
    Eigen::MatrixXd const & computeNonlinearMatrix();
    Eigen::MatrixXd const & computeTransformMatrix();    
    Eigen::MatrixXd const & computePsfMatrix();

    std::pair<int, int> const getImageDimensions() const {
        return _imageDimensions;
    }
    int const & getImageWidth() const {return _imageDimensions.second;}
    int const & getImageHeight() const {return _imageDimensions.first;}
    int const getImageSize() const {
        return getImageHeight()*getImageWidth();
    }
    virtual int getLinearSize() const = 0;
    virtual int getNonlinearSize() const = 0;
    virtual int getPsfBasisSize() const {return 0;}
    
    virtual Coordinate getCenter() const = 0;
    
    virtual ~Model(){}
protected:
    virtual void updateParametrizedImage() = 0;
    virtual void updateConstantImage() {}
    virtual void updateLinearMatrix()=0;    
    virtual void updateNonlinearMatrix()=0;
    virtual void updateTransformMatrix() =0;    
    virtual void updatePsfMatrix() {}

    std::pair<int,int> imageDimensions;

    VectorXd _parameterizedImage;
    VectorXd _constantImage;
    MatrixXd _linearMatrix;
    MatrixXd _nonlinearMatrix;
    MatrixXd _transformMatrix;
    MatrixXd _psfMatrix;

    Model(Model::Ptr const & other) {}

    Model(int const imageHeight,
            int const imageWidth,
            int const nonlinearSize,
            int const linearSize,
            int const psfBasisSize = 0,
            int const transformSize = 6)
        : _imageDimensions(imageHeight, imageWidth) {
        init(nonlinearSize, linearSize, psfBasisSize, transformSize);        
    }
    
private:
    void operator=(Model const & other) {}

    void init(int const nonlinearSize,
            int const linearSize,
            int const psfBasisSize,
            int const transformSize);
};


}}} //end namespace lsst::meas::multifit
#endif
