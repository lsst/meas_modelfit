#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include "boost/shared_ptr.hpp"

#include "Eigen/Core"
#include "ndarray/ndarray.hpp"

namespace lsst {
namespace meas {
namespace multifit {

typedef Eigen::Vector2d Coordinate;
typedef ndarray::ArrayRef<double, 1, 1> ParameterVector;
//dimensions are pixel height, pixel width
typedef ndarray::ArrayRef<double, 2, 2> ImageVector;
// dimensions are parameters, pixel height, pixel width
typedef ndarray::ArrayRef<double, 3, 3> DerivativeMatrix;

inline Eigen::Map<Eigen::VectorXd> extractEigenView(ParameterVector const & array) {
    return array.core().vector();
}

inline Eigen::Map<Eigen::VectorXd> extractEigenView(ImageVector const & array) {
    return array.core().vector();
}

inline Eigen::Map<Eigen::MatrixXd> extractEigenView(DerivativeMatrix const & array) {
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
    VectorXd const & getLinearParameters() {return _linearParameters;}
    VectorXd const & getNonlinearParameters() {return _nonlinearParameters;}
    /**
     * Apply transform to the model "after" any existing transform.
     */
    virtual void addTransform(Eigen::Transform2d const & transform) = 0;
    
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

    std::pair<int, int> const getImageDimensions() const;
    int const getImageWidth() const;
    int const getImageHeight() const;

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

    VectorXd _linearParameters;
    VectorXd _nonlinearParameters;
    VectorXd _modelImage;
    VectorXd _constantImage;
    MatrixXd _linearMatrix;
    MatrixXd _nonlinearMatrix;
    MatrixXd _transformMatrix;
    MatrixXd _psfMatrix;

    Model(Model const & other)
    {}
    Model(int const imageWidth,
            int const imageHeight,
            int const linearSize,
            int const nonlinearSize,
            int const psfBasisSize
            int const transformSize = 6);
private:
    void operator=(Model const & other) {}
};


}}} //end namespace lsst::meas::multifit
#endif
