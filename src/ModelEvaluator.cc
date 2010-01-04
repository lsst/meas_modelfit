#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/afw/image/MaskedImage.h"

namespace multifit = lsst::meas::multifit;

ndarray::Array<multifit::Pixel const, 1, 1> multifit::ModelEvaluator::computeModelImage() {
    if(!(_validProducts & MODEL_IMAGE)) {
        ProjectionFrameIterator i(_projectionList.begin());
        ProjectionFrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->computeModelImage();
        }
    }    
    return _modelImage;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeLinearParameterDerivative() {
    if (!(_validProducts & LINEAR_PARAMETER_DERIVATIVE)) {
        ProjectionFrameIterator i(_projectionList.begin());
        ProjectionFrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->computeLinearParameterDerivative();
        }
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }    
    return _linearParameterDerivative;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeNonlinearParameterDerivative() {
    if(!(_validProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        ProjectionFrameIterator i(_projectionList.begin());
        ProjectionFrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->computeNonlinearParameterDerivative();
        }
        _validProducts |= NONLINEAR_PARAMETER_DERIVATIVE;
    }    
    return _nonlinearParameterDerivative;
}


