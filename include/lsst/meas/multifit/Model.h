#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include "boost/shared_ptr.hpp"

#include "Eigen/Core"
#include "ndarray/ndarray.hpp"

namespace lsst {
namespace meas {
namespace multifit {
typedef ndarray::ArrayRef<double, 1, 1> ParameterVector;
typedef ndarray::ArrayRef<double, 2, 2> ImageVector;
typedef ndarray::ArrayRef<double, 3, 2> DerivativeMatrix;

class TransformedModel;

class Model {
public:
    typedef boost::shared_ptr<Model> Ptr;
    typedef boost::shared_ptr<const Model> ConstPtr;

    virtual TransformedModel * transform(Eigen::Transform2d transform, Psf::Ptr psf) const;

    virtual void evalImage(ParameterVector const & linearParameters, 
            ParameterVector const & nonlinearParameters,
            ImageVector const & image           
    );
    virtual void evalLinearDerivative(ParameterVector const & nonlinearParameters,
            DerivativeMatrix const & linearDerivative,
            ImageVector const & constantModel
    ) = 0;
    virtual void evalNonlinearDerivative(
            ParameterVector const & linearParameters, 
            ParameterVector const & nonlinearParameters,
            DerivateMatrix const & nonlinearDerivative,
    ) = 0;


    virtual int getNumLinearParameters() const = 0;
    virtual int getNumNonlinearParameters() const = 0;
    
protected:

    Model(Model const & other)
    {}

    virtual void transformNonlinear(
            Eigen::Transform2d const & transform, 
            ParameterVector const & input,
            ParameterVector const & output
    );
    
    virtual void evalTransformDerivative_impl(
            Eigen::Transform2d const & transform,
            ParameterVector const & linearParameters, 
            ParameterVector const & nonlinearParameters,
            DerivateMatrix const & transformDerivative
    );

private:
    void operator=(Model const & other) {}
};


class TransformedModel : public Model {
public:
    TransformedModel(
            Model::ConstPtr model, 
            Eigen::Transform2d transform, 
            Psf::Ptr psf,
            int overSampling = 1
    ) : _model(model), _transform(transform), _psf(psf) {
        
    }


    virtual void evalPsfDerivative(
            ParameterVector const & linearParameters, 
            ParameterVector const & nonlinearParameters,
            DerivateMatrix const & psfDerivative
    );

    virtual void evalTransformDerivative(
            ParameterVector const & linearParameters, 
            ParameterVector const & nonlinearParameters,
            DerivateMatrix const & transformDerivative
    ) {
        _model->evalTransformDerivative_impl(_transform, linearParameters,
                nonlinearParameters, transformDerivative);        
    }

    int getNumPsfParameters() const {return _psf->getNumParams();}

protected:
    TransformModel(TransformModel const & other)
        : _model(other._model), _transform(other._transform), _psf(other._psf) 
    {}
    Model::ConstPtr _model;
    Eigen::Transform2d _transform;
    Psf::Ptr _psf;
}


}}} //end namespace lsst::meas::multifit
#endif
