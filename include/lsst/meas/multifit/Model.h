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


    // should add to existing image rather than overwrite it
    virtual void evalParametrizedImage(ImageVector const & parametrizedImage) const = 0 ;
    // should add to existing image rather than overwrite it
    // by default, there is nothing to do. Only constrained models will need to
    // override this method
    virtual void evalConstantImage(ImageVector const& constantImage) const {};
    // CAUTION: Assumes DerivativeMatrix is zeroed. Will overwrite
    virtual void evalLinearDerivative(DerivativeMatrix const & linearDerivative) const = 0;    
    // CAUTION: Assumes DerivativeMatrix is zeroed. Will overwrite
    virtual void evalNonlinearDerivative(DerivateMatrix const & nonlinearDerivative) const = 0;
    // CAUTION: Assumes DerivativeMatrix is zeroed. Will overwrite
    virtual void evalTransformDerivative(DerivativeMatrix const & transformDerivative) const = 0;    
    // CAUTION: Assumes DerivativeMatrix is zeroed. Will overwrite
    virtual void evalPsfDerivative(DerivativeMatrix const & psfDerivative) const = 0;    


    virtual int getNumLinearParameters() const = 0;
    virtual int getNumNonlinearParameters() const = 0;
    /**
     * affine transforms have exactly 6 parameters. yay for magic numbers
     */
    virtual int getNumTransformParameters() const {return 6;}
    virtual int getNumPsfParameters() const {return 0;}
    
    virtual Coordinate getCenter() const = 0;
    
    virtual ~Model(){}
protected:

    Model(Model const & other)
    {}

private:
    void operator=(Model const & other) {}
};


}}} //end namespace lsst::meas::multifit
#endif
