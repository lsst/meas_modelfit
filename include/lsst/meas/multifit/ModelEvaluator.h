// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H
#define LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H

#include <vector>
#include <iostream>
#include "Eigen/Core"
#include "Eigen/LU"

#include "ndarray.hpp"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/meas/multifit/CharacterizedExposure.h"

namespace lsst {
namespace meas {
namespace multifit{

class ModelEvaluator : private boost::noncopyable {
public:

#ifndef SWIG
    class ProjectionFrame {
    public:
        ProjectionFrame(){}
        ProjectionFrame(ProjectionFrame const & other) : _projection(other._projection) {}
        explicit ProjectionFrame(ModelProjection::Ptr const & projection) : _projection(projection) {} 

        ModelProjection::ConstPtr const getModelProjection() const {return _projection;}
        WcsConstPtr const & getWcs() const {return _projection->getWcs();}
        FootprintConstPtr const & getFootprint() const {return _projection->getFootprint();}
        
        ndarray::Array<Pixel const, 1, 1> const getImageVector() const {return _imageVector;}
        ndarray::Array<Pixel const, 1, 1> const getVarianceVector() const {return _varianceVector;}
        ndarray::Array<Pixel const, 1, 1> const computeModelImage() {
            return _projection->computeModelImage();
        }
        ndarray::Array<Pixel const, 2, 1> const computeLinearParameterDerivative() {
            return _projection->computeLinearParameterDerivative();
        }
        ndarray::Array<Pixel const, 2, 1> const computeNonlinearParameterDerivative() {
            return _projection->computeNonlinearParameterDerivative();
        }

    private:
        friend class ModelEvaluator;

        ModelProjection::Ptr _projection;
        ndarray::Array<Pixel, 1, 1> _imageVector;
        ndarray::Array<Pixel, 1, 1> _varianceVector;
    };
    typedef std::list<ProjectionFrame> ProjectionFrameList;

#endif
    
    explicit ModelEvaluator(Model::ConstPtr const & model, int const nMinPix=-1) 
      : _validProducts(0),
        _model(model->clone())
    {
        setMinPixels(nMinPix);
    }

    template <typename ImagePixel, typename MaskPixel, typename VariancePixel>
    ModelEvaluator(
        Model::ConstPtr const & model, 
        std::list< boost::shared_ptr< CharacterizedExposure<
                ImagePixel, MaskPixel, VariancePixel
            > > > const & exposureList,
        int const nMinPix = 0 
    ) : _validProducts(0),
        _model(model->clone()) 
    {        
        setMinPixels(nMinPix);
        setExposureList(exposureList);
    }
    template<typename ImagePixel, typename MaskPixel, typename VariancePixel>
    void setExposureList(
        std::list< boost::shared_ptr< CharacterizedExposure<
                ImagePixel, MaskPixel, VariancePixel
            > > > const & exposureList
    );

#ifndef SWIG
    ndarray::Array<Pixel const, 1, 1> getImageVector() const {return _imageVector;}
    ndarray::Array<Pixel const, 1, 1> getVarianceVector() const {return _varianceVector;}
    ndarray::Array<Pixel const, 1, 1> computeModelImage();
    ndarray::Array<Pixel const, 2, 2> computeLinearParameterDerivative();
    ndarray::Array<Pixel const, 2, 2> computeNonlinearParameterDerivative();
#endif

    int const & getMinPixels() const {return _nMinPix;}
    void setMinPixels(int const nMinPix) {
        _nMinPix = (nMinPix > 0)? nMinPix : 0;
    }

    int const getLinearParameterSize() const {
        return _model->getLinearParameterSize();
    }
    int const getNonlinearParameterSize() const {
        return _model->getNonlinearParameterSize();
    }
    ParameterVector const & getLinearParameters() const {
        return _model->getLinearParameters();
    }
    ParameterVector const & getNonlinearParameters() const {
        return _model->getNonlinearParameters();
    }

    void setLinearParameters(
        ParameterConstIterator const & parameterIterator
    ) {    
        _model->setLinearParameters(parameterIterator);
        _validProducts &= (~MODEL_IMAGE);
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }
    void setNonlinearParameters(
        ParameterConstIterator const & parameterIterator
    ) {
        _model->setNonlinearParameters(parameterIterator);
        _validProducts &= (~MODEL_IMAGE);
        _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }
   
    Model::ConstPtr getModel() const {return _model;}
    int const getNProjections() const {return _projectionList.size();}
    int const getNPixels() const {return _imageVector.getSize<0>();}

#ifndef SWIG
    ProjectionFrameList const & getProjectionList() const {
        return _projectionList;
    }
#endif

private:        
    typedef ProjectionFrameList::iterator ProjectionFrameIterator;
    typedef Footprint::SpanList SpanList;

    enum ProductFlag {
        MODEL_IMAGE = 1<<0,
        LINEAR_PARAMETER_DERIVATIVE = 1<<1,
        NONLINEAR_PARAMETER_DERIVATIVE = 1<<2,
    };

    int _nMinPix;
    int _validProducts;
    Model::Ptr _model;    
    ProjectionFrameList _projectionList;
    
    ndarray::Array<Pixel, 1, 1> _imageVector;
    ndarray::Array<Pixel, 1, 1> _varianceVector;
    ndarray::Array<Pixel, 1, 1> _modelImage;
    ndarray::Array<Pixel, 2, 2> _linearParameterDerivative;
    ndarray::Array<Pixel, 2, 2> _nonlinearParameterDerivative;
};

}}} //end namespace lsst::meas::multifit

#endif
