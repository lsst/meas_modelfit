#ifndef LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H
#define LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "ndarray_fwd.hpp"

#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/projections/ModelProjection.h"

namespace lsst {
namespace meas {
namespace multifit{

class ModelEvaluator {
public:  
    typedef boost::shared_ptr<projections::ModelProjection> ModelProjectionPtr;
    typedef boost::shared_ptr<Kernel const> KernelPtr;
    typedef boost::shared_ptr<Wcs const> WcsPtr;
    typedef boost::shared_ptr<Footprint const> FootprintPtr;
    typedef boost::shared_ptr<Exposure const> ExposurePtr;
    typedef boost::shared_ptr<MaskedImage const> MaskedImagePtr;

    typedef boost::tuple<ExposurePtr, KernelPtr, double> CalibratedExposure;
    typedef std::list<CalibratedExposure> CalibratedExposureList;

    struct ProjectionFrame {        
        Frame(){}
        Frame(Frame const & other) 
          : projection(other.projection),
            exposure(other.exposure),
            kernel(other.kernel),
            wcs(other.wcs),
            footprint(other.footprint)
        {}
        Frame(
            ModelProjectionPtr const & projectionPtr, 
            ExposurePtr const & exposurePtr, 
            KernelPtr const & kernelPtr, 
            WcsPtr const & wcsPtr, 
            FootprintPtr const &footprintPtr
        ) : projection(projectionPtr), 
            exposure(exposurePtr),
            kernel(kernelPtr),
            wcs(wcsPtr),
            footprint(footprintPtr)
        {}

        void compressExposure();

        ModelProjectionPtr const projection;
        ExposurePtr const exposure;
        KernelPtr const kernel;
        WcsPtr const wcs;
        FootprintPtr const footprint;

        //references into larger matrices
        ndarray::Array<const Pixel, 1, 1> const imageVector;
        ndarray::Array<const Pixel, 1, 1> const varianceVector;
        ndarray::Array<Pixel, 1, 1> const modelImage;
        ndarray::Array<Pixel, 2, 2> const linearParameterDerivative;
        ndarray::Array<Pixel, 2, 2> const nonlinearParameterDerivative;
        ndarray::Array<Pixel, 2, 2> const wcsParameterDerivative;
        ndarray::Array<Pixel, 2, 2> const psfParameterDerivative;
    };
     
    typedef std::list<ProjectionFrame> ProjectionFrameList;

    explicit ModelEvaluator(Model::ConstPtr model, int activeProducts=0) :    
        _activeProducts(activeProducts), 
        _validProducts(0), 
        _model(model->clone())
    {}

    ModelEvaluator(
        Model::ConstPtr model, 
        CalibratedExposureList const & exposureList, 
        int activeProducts = 0
    ) : _activeProducts(activeProducts), 
        _validProducts(0), 
        _model(model->clone()) 
    {
        setExposureList(exposureList);
    }

    void setExposureList(
        CalibratedExposureList const & exposureList, 
        int nMinPix = 0
    );

    ndarray::Array<Pixel const, 1, 1> computeModelImage();
    ndarray::Array<Pixel const, 2, 2> computeLinearParameterDerivative();
    ndarray::Array<Pixel const, 2, 2> computeNonlinearParameterDerivative();
    ndarray::Array<Pixel const, 2, 2> computeWcsParameterDerivative();
    ndarray::Array<Pixel const, 2, 2> computePsfParameterDerivative();

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
    }
    void setNonlinearParameters(
        ParameterCosntIterator const & parameterIterator
    ) {
        _model->setNonlinearParameters(parameterIterator);
    }
    
    Model::ConstPtr getModel() const {return _model;}
    int const getNProjections() const {return _projectionList.size();}

    FrameList const & getProjectionList() const {return _projectionList;}

private:        
    typedef ProjectionFrameList::iterator ProjectionFrameIterator;    
    typedef CalibratedExposureList::const_iterator CalibratedExposureIterator;
    typedef Footprint::SpanList SpanList;

    int _activeProducts;
    int _validProducts;
    Model::ConstPtr _model;    

    ProjectionFrameList _projectionList;

    static FootprintPtr fixFootprint(FootprintPtr const &, MaskedImagePtr const &);

    ndarray::Array<Pixel, 1, 1> _imageVector;
    ndarray::Array<Pixel, 1, 1> _varianceVector;
    ndarray::Array<Pixel, 1, 1> _modelImage;
    ndarray::Array<Pixel, 2, 2> _linearParameterDerivative;
    ndarray::Array<Pixel, 2, 2> _nonlinearParameterDerivative;
    ndarray::Array<Pixel, 2, 2> _wcsParameterDerivative;
    ndarray::Array<Pixel, 2, 2> _psfParameterDerivative;
};

}}} //end namespace lsst::meas::multifit

#endif
