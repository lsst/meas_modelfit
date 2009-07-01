#ifndef LSST_MEAS_MULTIFIT_GAUSSIAN_PROBABILITY_H
#define LSST_MEAS_MULITFIT_GAUSSIAN_PROBABILITY_H

#include <set>

#include "Eigen/Core"
#include "lsst/meas/multifit/Model.h"

namespace lsst {
namespace meas {
namespace multifit {

struct ProbabilityExpansion {
    double _p;
    Eigen::VectorXd _dp;
    Eigen::MatrixXd _ddp;
    
    template <typename VecT, typename MatT>
    ProbabilityExpansion(double p, Eigen::MatrixBase<VecT> const & dp,
                         Eigen::MatrixBase<MatT> const& ddp)
        : _p(p), _dp(dp), _ddp(ddp) {}
};

class Probability {
    ObjectModel::Map _model_map;
public:
    virtual Eigen::VectorXd max() const = 0;
    virtual ProbabilityExpansion at(Eigen::VectorXd const & point) const = 0;
    const ObjectModel::Map& getModels() const { return _model_map; }
    
    // marginalizes over models not contained in tags
    virtual Probability* extract(std::set<ObjectModel::Key> const & tags) const = 0;

    // default implementation constructs tags set and calls the above
    virtual Probability* extract(ObjectModel::Key objectId) const;
    virtual ~Probability(){}
};

class GaussianProbability : public Probability {
    ProbabilityExpansion _expansion;
    Eigen::MatrixXd _covariance;
public:

    virtual Eigen::VectorXd max() const { return _expansion._dp; }

    virtual ProbabilityExpansion at(Eigen::VectorXd const & point) const { 
        return _expansion; 
    }
    
    // marginalizes over models not contained in tags
    virtual GaussianProbability* extract(std::set<ObjectModel::Key> const & tags) const;
    
    virtual GaussianProbability* extract(ObjectModel::Key objectId) const {
        return static_cast<GaussianProbability*>(
                Probability::extract(objectId));
    }
    
    Eigen::MatrixXd const & getCovariance() const {return _covariance;}
};

}}} //end namespace lsst::meas::multifit
#endif

