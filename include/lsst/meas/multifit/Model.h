#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <Eigen/Core>

#include <ndarray_fwd.hpp>

#include <lsst/afw/math/AffineTransform.h>
#include <lsst/afw/image/Wcs.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/detection/Footprint.h>

namespace lsst {
namespace meas {
namespace multifit {

class Model {
public:
    typedef lsst::afw::math::Kernel Kernel;
    typedef lsst::afw::detection::Footprint Footprint;
    typedef boost::shared_ptr<const Footprint> FootprintConstPtr;
    typedef lsst::afw::image::Wcs Wcs;
    typedef boost::shared_ptr<const Wcs> WcsConstPtr;


    typedef boost::shared_ptr<Model> Ptr;
    typedef boost::shared_ptr<const Model> ConstPtr;

    typedef double Pixel;
    typedef double * ParameterIterator;
    typedef double const * ParameterConstIterator;

    class Factory;
    class Definition;

    enum ProductFlag {
        NONE                             = 0x00,
        MODEL_IMAGE                      = 0x01,
        LINEAR_PARAMETER_DERIVATIVE      = 0x02,
        NONLINEAR_PARAMETER_DERIVATIVE   = 0x04,
        OBJECT_PARAMETER_DERIVATIVE      = 0x06,
        WCS_PARAMETER_DERIVATIVE         = 0x08,
        PSF_PARAMETER_DERIVATIVE         = 0x10,
        CALIBRATION_PARAMETER_DERIVATIVE = 0x18,
        ALL                              = 0x1F
    };

    virtual ~Model() {}

    virtual Model * clone() const = 0;

    /**
     *  Return a new Model of the same type and same parameters with a new coordinate system
     *  and Kernel.
     */
    virtual void reproject(
        Kernel const & kernel,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint,
        double photFactor
    ) = 0;

    Eigen::VectorXd const & getLinearParameters() const;
    void getLinearParameters(ParameterIterator parameters) const;
    void setLinearParameters(ParameterConstIterator parameters);

    Eigen::VectorXd const & getNonlinearParameters() const;
    void getNonlinearParameters(ParameterIterator parameters) const;
    virtual void setNonlinearParameters(ParameterConstIterator parameters);

    int const getLinearParameterSize() const;
    int const getNonlinearParameterSize() const;
    virtual int const getWcsParameterSize() const = 0;
    virtual int const getPsfParameterSize() const = 0;

    /**
     *  Compute the model image, flattened via its footprint.
     */
    virtual void computeModelImage(ndarray::Array<Pixel,1,1> const & vector) = 0;

    /**
     *  Compute the derivative of the model image with respect to its linear 
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computeLinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    ) = 0;

    /**
     *  Compute the derivative of the model image with respect to its nonlinear
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computeNonlinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    ) = 0;

    /**
     *  Compute the derivative of the model image with respect to its Wcs 
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computeWcsParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    ) = 0;

    /**
     *  Compute the derivative of the model image with respect to its PSF 
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computePsfParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    ) = 0;

    void enableProducts(int toAdd) {
        _activeProducts |= _enableProducts(toAdd);
    }

    void disableProducts(int toRemove) {
        _activeProducts &= (~_disableProducts(toRemove));
    }

    int const getActiveProducts() const { return _activeProducts; }
    double getPhotFactor() const { return _photFactor; }
    WcsConstPtr getWcs() const { return _wcs; }
    FootprintConstPtr getFootprint() const { return _footprint; }

    Definition const & getDefinition() const;

    boost::shared_ptr<Factory const> getFactory() const;

protected:

    Model(
        Definition const & definition,
        int activeProducts,
        WcsConstPtr const &wcs,
        FootprintConstPtr const & footprint,
        double photFactor
    );

    Model(Model const & other);

    void setProjectionVariables(
        WcsConstPtr wcs,
        FootprintConstPtr const & footprint,
        double photFactor
    ) {
        _wcs = wcs;
        _footprint = footprint;
        _photFactor = photFactor;
    }

    virtual int _enableProducts(int toAdd) = 0;
    virtual int _disableProducts(int toRemove) = 0;

    virtual void _handleLinearParameterChange() = 0;
    virtual void _handleNonlinearParameterChange() = 0;

private:

    int _activeProducts;
    double _photFactor;
    boost::scoped_ptr<Definition> _definition;
    FootprintConstPtr _footprint;
    WcsConstPtr _wcs;

    void operator=(Model const & other) {} // assignment is disabled

};

class Model::Factory {
public:
    typedef boost::shared_ptr<Factory> Ptr;
    typedef boost::shared_ptr<const Factory> ConstPtr;
    
    virtual Model * project(
        ParameterConstIterator linearParameterBegin,
        ParameterConstIterator const linearParameterEnd,
        ParameterConstIterator nonlinearParameterBegin,
        Kernel const & kernel,
        WcsConstPtr const &wcs,
        FootprintConstPtr const &footprint,
        double photFactor
    ) const = 0;

    virtual int const getNonlinearParameterSize() const = 0;

    virtual int const getMinLinearParameterSize() const = 0;
    virtual int const getMaxLinearParameterSize() const = 0;

    virtual ~Factory() {}
};

class Model::Definition {
public:
    typedef boost::shared_ptr<Factory> Ptr;
    typedef boost::shared_ptr<const Factory> ConstPtr;
    
    Eigen::VectorXd const & getLinearParameters() const { 
        return _linearParameters; 
    }
    Eigen::VectorXd const & getNonlinearParameters() const { 
        return _nonlinearParameters; 
    }
    
    Factory::ConstPtr const & getFactory() const { return _factory; }

    void setLinearParameters(ParameterConstIterator parameters) {
        std::copy(
                parameters, 
                parameters+getLinearParameterSize(), 
                _linearParameters.data()
        );
    }

    void setNonlinearParameters(ParameterConstIterator parameters) {
        std::copy(
                parameters, 
                parameters+getNonlinearParameterSize(), 
                _nonlinearParameters.data()
        );
    }

    int const getLinearParameterSize() const { 
        return _linearParameters.size(); 
    }
    int const getNonlinearParameterSize() const { 
        return _nonlinearParameters.size(); 
    }
    
    Model * project(
            Kernel const & kernel,
            WcsConstPtr const &wcs,
            FootprintConstPtr const &footprint,
            double photFactor
    ) const {
        return _factory->project(
                _linearParameters.data(),
                _linearParameters.data() + getLinearParameterSize(),
                _nonlinearParameters.data(),
                kernel, wcs, footprint, photFactor
        );
    }

    explicit Definition(
            Factory::ConstPtr const & factory, 
            int linearParameterSize
    ) : _factory(factory),
        _linearParameters(linearParameterSize),
        _nonlinearParameters(factory->getNonlinearParameterSize())
    {}
    
    Definition(
        Factory::ConstPtr const & factory,
        ParameterConstIterator linearParameterBegin,
        ParameterConstIterator const linearParameterEnd,
        ParameterConstIterator nonlinearParameterBegin
    );

private:
    Factory::ConstPtr _factory;
    Eigen::VectorXd _linearParameters;
    Eigen::VectorXd _nonlinearParameters;
};

inline Eigen::VectorXd const & Model::getLinearParameters() const {
    return _definition->getLinearParameters(); 
}

inline void Model::setLinearParameters(ParameterConstIterator parameters) {
    _definition->setLinearParameters(parameters);
    _handleLinearParameterChange();
}

inline Eigen::VectorXd const & Model::getNonlinearParameters() const {
    return _definition->getNonlinearParameters(); 
}

inline void Model::setNonlinearParameters(ParameterConstIterator parameters) {
    _definition->setNonlinearParameters(parameters);
    _handleNonlinearParameterChange();
}

inline int const Model::getLinearParameterSize() const {
    return (_activeProducts & LINEAR_PARAMETER_DERIVATIVE) ? getLinearParameters().size() : 0; 
}

inline int const Model::getNonlinearParameterSize() const {
    return (_activeProducts & NONLINEAR_PARAMETER_DERIVATIVE) ? getNonlinearParameters().size() : 0;
}

inline Model::Definition const & Model::getDefinition() const { return *_definition; }

inline Model::Factory::ConstPtr Model::getFactory() const { return _definition->getFactory(); }

}}} //namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_H
