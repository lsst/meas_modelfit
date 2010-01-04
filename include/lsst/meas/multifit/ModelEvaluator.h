#ifndef LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H
#define LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "ndarray.hpp"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/footprintUtils.h"

namespace lsst {
namespace meas {
namespace multifit{

class ModelEvaluator {
public:     
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
   
    template<typename ImagePixel, typename MaskPixel, typename VariancePixel> 
    class Traits{
    public:                
        typedef typename lsst::afw::image::Exposure<ImagePixel, MaskPixel, VariancePixel> Exposure;
        typedef typename boost::shared_ptr<Exposure const> ExposureConstPtr;
        typedef typename std::pair<ExposureConstPtr, PsfConstPtr> CalibratedExposure;
        typedef typename std::list<CalibratedExposure> CalibratedExposureList;
        typedef typename CalibratedExposureList::const_iterator CalibratedExposureIterator;
        typedef typename std::list<ExposureConstPtr> ExposureList;

    };        

    typedef std::list<ProjectionFrame> ProjectionFrameList;

    explicit ModelEvaluator(Model::ConstPtr const & model, int const nMinPix=-1) 
      : _validProducts(0),
        _model(model->clone())
    {
        setMinPixels(nMinPix);
    }

    template<typename ImagePixel, typename MaskPixel, typename VariancePixel>
    ModelEvaluator(
        Model::ConstPtr const & model, 
        typename Traits<ImagePixel, MaskPixel, VariancePixel>::CalibratedExposureList const & exposureList,
        int const nMinPix = 0 
    ) : _validProducts(0),
        _model(model->clone()) 
    {        
        setMinPixels(nMinPix);
        setExposureList<ImagePixel, MaskPixel, VariancePixel>(exposureList);
    }

    int const & getMinPixels() const {return _nMinPix;}
    void setMinPixels(int const nMinPix) {
        _nMinPix = (nMinPix > 0)? nMinPix : 0;
    }

    template<typename ImagePixel, typename MaskPixel, typename VariancePixel>
    void setExposureList(
        typename Traits<ImagePixel, MaskPixel, VariancePixel>::CalibratedExposureList const & exposureList
    ) {
        typedef Traits<ImagePixel, MaskPixel, VariancePixel> Traits;
        
        _projectionList.clear();
        _validProducts = 0;

        int nLinear = getLinearParameterSize();
        int nNonlinear = getNonlinearParameterSize();

        int pixSum = 0;

        ModelProjection::Ptr projection;
        typename Traits::ExposureConstPtr exposure;
        FootprintConstPtr footprint;
        PsfConstPtr psf;
        WcsConstPtr wcs;    
  
        //exposures which contain fewer than _nMinPix pixels will be rejected
        //construct a list containing only those exposure which were not rejected
        typename Traits::ExposureList goodExposureList;

        // loop to create projections
        for(typename Traits::CalibratedExposureIterator i(exposureList.begin()), end(exposureList.end());
            i != end; ++i
        ) {
            exposure = i->first;
            psf = i->second;
            wcs = exposure->getWcs();        
            footprint = _model->computeProjectionFootprint(psf, wcs);
            footprint = clipAndMaskFootprint<MaskPixel>(
                footprint, exposure->getMaskedImage().getMask()
            );

            //ignore exposures with too few contributing pixels        
            if (footprint->getNpix() > _nMinPix) {
                ProjectionFrame frame(
                    _model->makeProjection(psf, wcs, footprint)
                );
      
                _projectionList.push_back(frame);
                goodExposureList.push_back(exposure);

                pixSum += footprint->getNpix();
            }
        }

        //  allocate matrix buffers
        ndarray::shallow(_imageVector) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
        ndarray::shallow(_varianceVector) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));

        ndarray::shallow(_modelImage) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
        ndarray::shallow(_linearParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(nLinear, pixSum)
        );
        ndarray::shallow(_nonlinearParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(nNonlinear, pixSum)
        );    
    
        int nPix;
        int pixelStart = 0, pixelEnd;

        typename Traits::ExposureList::const_iterator exposureIter(goodExposureList.begin());

        //loop to assign matrix buffers to each projection Frame
        for(ProjectionFrameList::iterator i(_projectionList.begin()), end(_projectionList.end()); 
            i != end; ++end
        ) {
            ProjectionFrame & frame(*i);
            nPix = frame.getFootprint()->getNpix();
            pixelEnd = pixelStart + nPix;

            // set image/variance buffers
            ndarray::shallow(frame._imageVector) = _imageVector[
                ndarray::view(pixelStart, pixelEnd)
            ];
            ndarray::shallow(frame._varianceVector) = _varianceVector[
                ndarray::view(pixelStart, pixelEnd)
            ];
        
            CompressFunctor<ImagePixel, MaskPixel, VariancePixel> functor(
                (*exposureIter)->getMaskedImage(), 
                frame._imageVector, 
                frame._varianceVector
            );
            functor.apply(*frame.getFootprint());    

        
            //set modelImage buffer
            frame._projection->setModelImageBuffer(
                _modelImage[ndarray::view(pixelStart, pixelEnd)]
            );
        
            //set linear buffer
            frame._projection->setLinearParameterDerivativeBuffer(
                _linearParameterDerivative[ndarray::view()(pixelStart, pixelEnd)]
            );
            //set nonlinear buffer
            frame._projection->setNonlinearParameterDerivativeBuffer(
                _nonlinearParameterDerivative[ndarray::view()(pixelStart, pixelEnd)]
            );

            pixelStart = pixelEnd;
            ++exposureIter;
        }   
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
    
    ndarray::Array<Pixel, 1, 1> _imageVector;
    ndarray::Array<Pixel, 1, 1> _varianceVector;
    ndarray::Array<Pixel, 1, 1> _modelImage;
    ndarray::Array<Pixel, 2, 2> _linearParameterDerivative;
    ndarray::Array<Pixel, 2, 2> _nonlinearParameterDerivative;
};

}}} //end namespace lsst::meas::multifit

#endif
