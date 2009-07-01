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
    ModelEvaluator(ObjectModel::Map const & model_set, 
                   ExposureContainer exposures,
                   int marginalizationFlags,
                   Probability const * prior);


    template <typename ExposureContainer>
    ModelEvaluator(ObjectModel::Map const & model_set,
                   ndarray::ArrayRef<double,1,1> const & linear,
                   ndarray::ArrayRef<double,1,1> const & nonlinear,
                   ExposureContainer exposures,
                   int marginalizationFlags,
                   Probability const * prior = NULL);

    ObjectModel::Map const& getModels() const {return _models;}

    // evaluate the likelihood of the model at the current values 
    //(only linear parameters are updated)
    ProbabilityExpansion const & evaluate();

    Eigen::VectorXd step();

    std::pair<double,GaussianProbability> run(int end_condition);

    GaussianProbability getProbability() const;
  
    ObjectModel& editModel(const ObjectModel::Key& tag);
    ObjectModel::Ptr const getModel(ObjectModel::Key const& id) const;

    int getNumTotalLinear() const {return _linear.size();}
    int getNumTotalNonLinear() const {return _nonlinear.size();}
    int getNumTotalPixels() const {return _numTotalPixels;}
    int getNumExposures() const {return _exposures.size();}
    int getNumObjects() const {return _models.size();}
    int getNumSources() const {return _workspace.size();}

    class SourceView {
        Eigen::Block<Eigen::MatrixXd> _gm;  // block of grand matrix
        Eigen::Block<Eigen::MatrixXd> _gv;  // block of grand data vector
        Eigen::MatrixXd _gr; // copy of residuals with this source un-subtracted

        // more data to come
        SourceView(const ObjectModel::Key& tag, 
                const CalibratedExposure::Ptr & exposure, 
                const std::string & constraint);
    
        friend class ModelEvaluator;

    public:
        std::pair<double,GaussianProbability> evaluate() const;
    };

private:    
    typedef std::map<CalibratedExposure::Ptr, int>::const_iterator 
            ExposureIndexIterator;
    typedef ObjectModel::SourceModel SourceModel;
    typedef std::vector<SourceModel*> SourceModelList;
    typedef std::vector<SourceModel*>::iterator SourceModelIterator;

    template <typename ExposureIterator>
    void setExposures(
            ExposureIterator const & begin, 
            ExposureIterator const & end);

    //begin stuff that wraps Exposure interface -------------------------------

    void buildBkgSubtracted();

    int getNumBkgParam(CalibratedExposure::Ptr const & exposure) const {
        // TODO: query exposure for background parameter size
        return testFlags(BKG) * (0);
    }

    int getNumPsfParam(CalibratedExposure::Ptr const & exposure) const {
        // TODO: query exposure for PSF parameter size
        return testFlags(PSF) * (0);
    }

    void computeBkgMatrix(CalibratedExposure::Ptr const & exposure) {
        int height = exposure->getHeight();
        int width = exposure->getWidth();
        ndarray::ArrayRef<double,3,2> ref(_calibration_matrix.data(),
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

    int _numTotalPixels;
    std::map<CalibratedExposure::Ptr, int> _exposures;
    ObjectModel::Map _models;
    SourceModelList _workspace;

    Eigen::VectorXd _linear, _nonlinear;
    Eigen::VectorXd _bkg_subtracted;
    Eigen::VectorXd _residuals;
    ProbabilityExpansion _posterior;

    Eigen::MatrixXd _linear_matrix, _nonlinear_matrix, _calibration_matrix;

    int _marginalizationFlags;

    inline bool testFlags(MarginalizationFlagsEnum f) const { 
        return _marginalizationFlags & f; 
    }

    // returned array is (linear)
    ndarray::ArrayRef<double,1,1> getLinearParameterSection(
            ObjectModelIterator & iter) {
        return ndarray::ArrayRef<double,1,1>(
                _linear.data() + iter.getLinearOffset(),
                (*iter)->getNumLinearParam());
    }

    // returned array is (nonlinear)
    ndarray::ArrayRef<double,1,1> getNonlinearParameterSection(
            ObjectModelIterator const & iter) {
        ObjectModel::Ptr model = *iter;
        return ndarray::ArrayRef<double,1,1>(
                _nonlinear.data() + iter.getNonlinearOffset(),
                model->getNumNonlinearParam());
    }
        
    // returned array is (y,x)
    ndarray::ArrayRef<double,2,2> getResidualsSection(
            ExposureIndexIterator const & exposure) {
        int height = exposure->first->getHeight();
        int width = exposure->first->getWidth();
        return ndarray::ArrayRef<double,2,2>(
                _residuals.data()+exposure->second,
                ndarray::make_index(height, width));
    }

    // returned array is (linear,y,x)
    ndarray::ArrayRef<double,3,2>
    getLinearMatrixSection(ObjectModelIterator const & modelIter, 
            ExposureIndexIterator const & exposureIter) {
        CalibratedExposure::Ptr exposure = exposureIter->first;
        int height = exposure->getHeight();
        int width = exposure->getWidth();
        int offset = modelIter.getLinearOffset()*_numTotalPixels 
                + exposureIter->second;
        ObjectModel::Ptr model = *modelIter;
        return ndarray::ArrayRef<double,3,2>(
            _linear_matrix.data() + offset,
            ndarray::make_index(model->getNumLinearParam(), height, width),
            ndarray::make_index(_numTotalPixels, width,1));
    }

    // returned array is (nonlinear,y,x)
    ndarray::ArrayRef<double,3,2> getNonlinearMatrixSection(
            ObjectModelIterator const & modelIter, 
            ExposureIndexIterator const & exposureIter
    ) {
        CalibratedExposure::Ptr exposure = exposureIter->first;
        int height = exposure->getHeight();
        int width = exposure->getWidth();
        int offset = modelIter.getNonlinearOffset()*_numTotalPixels
                + exposureIter->second;
        ObjectModel::Ptr model = *modelIter;

        return ndarray::ArrayRef<double,3,2>(
            _nonlinear_matrix.data() + offset,
            ndarray::make_index(model->getNumNonlinearParam(), height,width),
            ndarray::make_index(_numTotalPixels, width,1));
    }

    ndarray::ArrayRef <double, 3, 2> getPsfMatrix(
            ExposureIndexIterator const & exposureIter
    ) {
        CalibratedExposure::Ptr exposure = exposureIter->first;
        int height = exposure->getHeight();
        int width = exposure->getWidth();
        int offset = testFlags(BKG) * exposure->getNumBackgroundParam();
        int size = testFlags(PSF) * exposure->getNumPsfParam();
        return ndarray::ArrayRef<double,3,2>(
                _calibration_matrix.data() + offset*width*height,
                ndarray::make_index(size,height,width),
                ndarray::make_index(height*width,width,1));        
    }

};

}}} //end namespace lsst::meas::multifit

#endif
