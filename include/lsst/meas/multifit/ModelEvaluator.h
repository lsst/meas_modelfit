#ifndef LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H
#define LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "ndarray.hpp"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"

namespace lsst {
namespace meas {
namespace multifit{

class ModelEvaluator {
public:     
    class ProjectionFrame {        
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
        template<typename ImageT> friend class Traits;
        friend class ModelEvaluator;

        ModelProjection::Ptr _projection;
        ndarray::Array<Pixel, 1, 1> _imageVector;
        ndarray::Array<Pixel, 1, 1> _varianceVector;
    };
   
    template<typename ImageT> 
    class Traits{
        typedef typename lsst::afw::image::MaskedImage<ImageT> MaskedImage;
        typedef typename lsst::afw::image::Exposure<ImageT> Exposure;
        typedef typename boost::shared_ptr<Exposure> ExposureConstPtr;
        typedef typename boost::tuple<ExposureConstPtr, PsfConstPtr> CalibratedExposure;
        typedef typename std::list<CalibratedExposure> CalibratedExposureList;
        
        static void setExposureList(ModelEvaluator &, CalibratedExposureList const &);
        static void compressExposure(ProjectionFrame &, ExposureConstPtr const &);
        static FootprintConstPtr fixFootprint(
            FootprintConstPtr const &, typename MaskedImage::MaskPtr const &
        );
    };

    typedef std::list<ProjectionFrame> ProjectionFrameList;

    explicit ModelEvaluator(Model::ConstPtr const & model, int const nMinPix=-1) 
      : _validProducts(0),
        _model(model->clone())
    {
        setMinPixels(nMinPix);
    }

    template<typename ImageT>
    ModelEvaluator(
        Model::ConstPtr const & model, 
        typename Traits<ImageT>::CalibratedExposureList const & exposureList,
        int const nMinPix = -1
    ) : _validProducts(0),
        _model(model->clone()) 
    {        
        setMinPixels(nMinPix);
        setExposureList<ImageT>(exposureList);
    }

    int const & getMinPixels() const {return _nMinPix;}
    void setMinPixels(int const nMinPix) {
        _nMinPix = (nMinPix > 0)? nMinPix : -1;
    }

    template <typename ImageT>
    void setExposureList(
        typename Traits<ImageT>::CalibratedExposureList const & exposureList
    ) {
        Traits<ImageT>::setExposureList(*this, exposureList);
    }
   
    ndarray::Array<Pixel const, 1, 1> getImageVector() const {return _imageVector;}
    ndarray::Array<Pixel const, 1, 1> getVarianceVector() const {return _varianceVector;}
    ndarray::Array<Pixel const, 1, 1> computeModelImage();
    ndarray::Array<Pixel const, 2, 2> computeLinearParameterDerivative();
    ndarray::Array<Pixel const, 2, 2> computeNonlinearParameterDerivative();


    int const getLinearParameterSize() const {
        return _model->getLinearParameterSize();
    }
    int const getNonlinearParameterSize() const {
        return _model->getNonlinearParameterSize();
    }
    ParameterVector const & getLinearParameters() const {
        return _model->getLinearParameterVector();
    }
    ParameterVector const & getNonlinearParameters() const {
        return _model->getNonlinearParameterVector();
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
    ProjectionFrameList const & getProjectionList() const {
        return _projectionList;
    }

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
    
        
    template<typename ImageT>
    static FootprintConstPtr fixFootprint(
        FootprintConstPtr const & footprint, 
        typename Traits<ImageT>::MaskedImage::MaskPtr const & mask
    ) {
        return Traits<ImageT>::fixFootprint(footprint, mask);
    }

    template<typename ImageT>
    static void compressExposure(
        ProjectionFrame & frame,
        typename Traits<ImageT>::ExposureConstPtr const & exposure
    ) {
        Traits<ImageT>::compressExposure(frame, exposure);
    }

    ndarray::Array<Pixel, 1, 1> _imageVector;
    ndarray::Array<Pixel, 1, 1> _varianceVector;
    ndarray::Array<Pixel, 1, 1> _modelImage;
    ndarray::Array<Pixel, 2, 2> _linearParameterDerivative;
    ndarray::Array<Pixel, 2, 2> _nonlinearParameterDerivative;
};

}}} //end namespace lsst::meas::multifit

#endif
