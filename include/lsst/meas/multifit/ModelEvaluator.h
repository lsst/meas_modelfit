#ifndef LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H
#define LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "ndarray/ndarray.hpp"

#include "GaussianProbability.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelIterator.h"
#include "lsst/meas/multifit/ConstrainedModel.h"
#include "lsst/meas/multifit/CalibratedExposure.h"

namespace lsst {
namespace meas {
namespace multifit{

class ModelEvaluator {
public:

    enum MarginalizationFlagsEnum { NONE=0, PSF=1, BKG=2 };

    template <typename ExposureContainer>
    ModelEvaluator(Model::ConstPtr model,
                   ExposureContainer exposures,
                   int marginalizationFlags,
                   Probability const * prior);


    template <typename ExposureContainer>
    ModelEvaluator(Model::ConstPtr model,
                   ndarray::ArrayRef<double,1,1> const & linear,
                   ndarray::ArrayRef<double,1,1> const & nonlinear,
                   ExposureContainer exposures,
                   int marginalizationFlags,
                   Probability const * prior = NULL);

    Model::ConstPtr getModel() const {return _model;}

    std::vector<TransformedModel::Ptr> const & getModelStack() const {return _modelStack;}

    // evaluate the likelihood of the model at the current values 
    //(only linear parameters are updated)
    ProbabilityExpansion const & evaluate();

    Eigen::VectorXd step();

    std::pair<double,GaussianProbability> run(int end_condition);

    GaussianProbability getProbability() const;
  
    void setModel(Model::ConstPtr model);
 
    int getNumLinear() const {return _linear.size();}
    int getNumNonlinear() const {return _nonlinear.size();}
    int getNumTotalPixels() const {return _numTotalPixels;}
    int getNumExposures() const {return _exposures.size();}

    class SourceView {
        Eigen::Block<Eigen::MatrixXd> _gm;  // block of grand matrix
        Eigen::Block<Eigen::MatrixXd> _gv;  // block of grand data vector
        Eigen::MatrixXd _gr; // copy of residuals with this source un-subtracted

        // more data to come
        SourceView(const CalibratedExposure::Ptr & exposure, 
                const std::string & constraint);
    
        friend class ModelEvaluator;

    public:
        std::pair<double,GaussianProbability> evaluate() const;
    };

private:    
    typedef std::vector<TransformedModel::Ptr> ModelStack;
    typedef ndarray::ArrayCore<double, 2, 2> ImageCore;
    typedef ndarray::ArrayCore<double, 3, 2> DerivativeCore;

    int _numTotalPixels;
    std::vector<CalibratedExposure::Ptr> _exposures;
    std::vector<ImageCore> _residualsSection;
    std::vector<DerivativeCore> _linearMatrixSection;
    std::vector<DerivativeCore> _nonlinearMatrixSection;

    Model::ConstPtr _model;    
    ModelStack _modelStack;


    Eigen::VectorXd _linear, _nonlinear;
    Eigen::VectorXd _bkgSubtracted;
    Eigen::VectorXd _residuals;
    ProbabilityExpansion _posterior;

    Eigen::MatrixXd _linearMatrix, _nonlinearMatrix, _calibrationMatrix;

    int _marginalizationFlags;

    template <typename ExposureContainer>
    void setExposures(ExposureContainer const & exposures);

    //begin stuff that wraps Exposure interface -------------------------------

    void buildBkgSubtracted();

    int getNumBkgParam(int const exposureId) const {
        return getNumBkgParam(_exposures.at(exposureId));
    }
    int getNumBkgParam(CalibratedExposure::Ptr const & exposure) const {
        return testFlags(BKG) * exposure->getNumBkgParam();
    }

    int getNumPsfParam(int const exposureId) const {
        return getNumPsfParam(_exposures.at(exposureId));
    }
    int getNumPsfParam(CalibratedExposure::Ptr const & exposure) const {
        return testFlags(PSF) * exposure->getNumPsfParam();
    }

    void computeBkgMatrix(CalibratedExposure::Ptr const & exposure) {
        int height = exposure->getHeight();
        int width = exposure->getWidth();
        ndarray::ArrayRef<double,3,2> ref(_calibrationMatrix.data(),
                ndarray::make_index(getNumBkgParam(exposure),height,width),
                ndarray::make_index(height*width,width,1));

        // TODO: query exposure for background matrix, append to ref
    }

    void addCalibrationFisherMatrix(CalibratedExposure::Ptr const & exposure, 
            Eigen::MatrixXd & matrix) const {
        // TODO: get exposure to add inverse calibration matrix to the 
        // given matrix.  Should only add blocks that correspond to the 
        // set bits in the marginalization flags.
    }

    //end stuff that wraps Exposure interface ---------------------------------



    inline bool testFlags(int flags) const { 
        return _marginalizationFlags & flag; 
    }
       
    // returned array is (y,x)
    ImageVector getResidualsSection(int exposureId) {
        return ndarray::ArrayRef<double, 2, 2>(residualsSection.at(exposureId);
    }

    // returned array is (linear,y,x)
    DerivateMatrix getLinearMatrixSection(int exposureId) {
        return DerivativeMatrix(linearMatrixSection.at(exposureId);
    }

    // returned array is (nonlinear,y,x)
    DerivativeMatrix getNonlinearMatrixSection(int exposureId) {
        return DerivativeMatrix(nonlinearMatrixSection.at(exposureId);
    }

    DerivativeMatrix getPsfMatrixSection(int exposureid) {
        return DerivativeMatrix(psfMatrixSection.at(exposureId);
    }

};

}}} //end namespace lsst::meas::multifit

#endif
