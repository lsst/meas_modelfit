#ifndef LSST_MEAS_MULTIFIT_CONVOLVEDMODEL_H
#define LSST_MEAS_MULTIFIT_CONVOLVEDMODEL_H

#include "ndarray/ndarray.hpp"
#include "Eigen/Core"
#include "lsst/meas/multifit/Model.h"

namespace lsst {
namespace meas {
namespace multifit {

class ConvolvedModel : public Model {
public:
    typedef boost::shared_ptr<ConvolvedModel> Ptr;
    typedef boost::shared_ptr<const ConvolvedModel> ConstPtr;

    ConvolvedModel(Model::ConstPtr unconvolved, Psf::ConstPtr psf, int oversample) :
            Model(unconvolved->getImageHeight(), unconvolved->getImageHeight(),
                    unconvolved->getNonlinearSize(), 
                    unconvolved->getLinearSize(), 
                    psf->getBasisSize()),
            _oversample((oversample > 0) ? oversample : 1),
            _unconvolved(unconvolved->project(
                    unconvolved->getImageHeight()*_oversample,
                    unconvolved->getImageWidth()*_oversample,
                    AffineTransform::makeScaling(_oversample)*unconvolved->getTransform())),
            _psf(psf),
    {
        setTransform(unconvolved->getTransform());
        _unconvolved->addTransform(AffineTransform::makeScaling(_oversample));
    }
    
    virtual ~ConvolvedModel() {}
    virtual Model * clone() const {
        return new ConvolvedModel(*this);
    }
    /**
     * Creates a convolved model
     * If the model already has a Psf, the old Psf is ignored.
     */
    virtual Model * convolve(Psf::ConstPtr psf) const {
        return new ConvolvedModel(_unconvolved, psf);
    };

    virtual void setLinearParameters(Eigen::VectorXd const & parameters) {
        _linearDirty = true;
        _constantDirty = true;
        _parametrizedDirty = true;
        _transformDirty = true;
        _psfDirty = true;
        _unconvolved->setLinearParameters(parameters);    
    }
    virtual void setNonlinearParameters(Eigen::VectorXd const & parameters) {
        _nonlinearDirty = true;
        _constantDirty = true;
        _parametrizedDirty = true;
        _transformDirty = true;
        _psfDirty = true;
        _unconvolved->SetNonlinearParameters(parameters);        
    }   
    Eigen::VectorXd void getNonLinearParameters() {
        return _unconvolved->getNonlinearParameters();    
    }
    Eigen::VectorXd void getLinearParameters() {
        return _unconvolved->getLinearParameters();
    }
    virtual void setTransform(AffineTransform const & transform) {
        _transformDirty = true;
        _transform = transform;    
        _unconvolved->setTransform(
                AffineTransform::makeScaling(_oversample)*transform
        );
    }
    virtual void addTransform(AffineTransform const & transform) {
        _transformDirty = true;
        _transform = transform * _transform;
        _unconvolved->addTransform(transform);        
    }
    
    virtual int getLinearSize() const {
        return _unconvolved->getLinearSize();
    }
    virtual int getNonlinearSize() const {
        return _unconvolved->getNonlinearSize();    
    }
    virtual int getPsfBasisSize() const {
        return _psfBasisSize();
    }
    
    virtual Coordinate getCenter() const = {return _unconvolved->center();}

    virtual void updateParametrizedImage();
    virtual void updateConstantImage();
    virtual void updateLinearMatrix();    
    virtual void updateNonlinearMatrix();
    virtual void updateTransformMatrix();    
    virtual void updatePsfMatrix();    
protected:
    explicit ConvolvedModel(ConvolvedModel const & other) 
        : _unconvolved(other._unconvolved->clone()), _psf(other._psf)
    {}

    int _oversample;
    Model::Ptr _unconvolved;
    Psf::ConstPtr _psf;   
    AffineTransform _transform;

    bool _constantDirty, _parametrizedDirty;
    bool _linearDirty, _nonlinearDirty;
    bool _transformDirty, _psfDirty;
};

}}} //namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_CONVOLVEDMODEL_H
