#ifndef LSST_MEAS_MULTIFIT_POINTSOURCEMODEL_H
#define LSST_MEAS_MULTIFIT_POINTSOURCEMODEL_H

#include "ndarray/ndarray.hpp"
#include "Eigen/Core"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/MultifitPsf.h"

namespace lsst {
namespace meas {
namespace multifit {

class PointSourceModel : public Model{
    typedef lsst::meas::multifit::MultifitPsf Psf;

    PointSourceModel(Coordinate const center, double linear, Psf::Ptr psf = null)
        : _center(nonlinear), _magnitude(linear), _psf(psf)
    {}
    ~PointSourceModel(){}

    virtual Model * clone() const {
        return new PointSourceModel(*this);
    }
    virtual Model * convolve(Psf::ConstPtr psf) {
        return new PointSourceModel(_center, _magnitude, psf);    
    }

    virtual void setNonlinearParameters(Eigen::VectorXd const & parameters) {
        if(parameters.size() != getNumLinearParameters())
            return;
        
        //deep copy, 2 elements
        _linear << parameters;
    }
    virtual void setLinearParameters(Eigen::VectorXd const & parameters) {
        if(parameters.size () <= getNumNonlinearParameters())
            return;

        //deep copy, 1 element
        _magnitude = parameters[0];
    }
    virtual void addTransform(Eigen::Transform2d const & transform) {
        _unconvolved->addTransform(transform);        
    }

    /**
     * PointSourceModel has exactly two nonlinear parameters x,y position    
     */
    virtual int getNumNonlinearParameters() const {return 2;}
    /** 
     * PointSourceModel has exactly one linear parameter
     */
    virtual int getNumLinearParamters() const {return 1;}
    virtual int getNumPsfParameters() const {
        return (_psf) ? _psf->getNumParam() : 0;
    }

    virtual Coordinate getCenter() {return _center;}

    virtual void evalParametrizedImage(ImageVector const & modelImage) const;
    virtual void evalLinearDerivative(DerivativeMatrix const & linearDerivative) const;
    virtual void evalNonlinearDerivative(DerivativeMatrix const & nonlinearDerivative) const;
    virtual void evalTransformDerivative(DerivativeMatrix const & transformDerivative) const;
    virtual void evalPsfDerivative(DerivativeMatrix const & psfDerivative) const;
    
protected:
    typedef Psf::ImageT ImageT;

    explicit PointSourceModel(PointSourceModel const & other) 
        : _nonlinear(other._nonlinear), _linear(other._linear), _psf(other._psf)
    {}

    Coordinate _center; //2 nonlinear parameters
    double _magnitude; //1 linear parameter

    Psf::Ptr _psf; 
    Eigen::Transform2d _transform;
};

}}} //namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_POINTSOURCEMODEL_H
