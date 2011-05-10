#include "lsst/meas/multifit/GaussNewtonOptimizer.h"
#include "lsst/meas/algorithms/shapelet/NLSolver.h"
#include "lsst/meas/multifit/Evaluation.h"

#include <Eigen/Core>
#include <cfloat>
#include "lsst/ndarray/eigen.h"

namespace {

class Solver : public lsst::meas::algorithms::shapelet::NLSolver {
public:
    Solver(lsst::meas::multifit::Evaluation & evaluation, 
           double const fTol, double const gTol, 
           double const minStep, int const maxIter, 
           double const tau
    ) : _nCoeff(evaluation.getEvaluator()->getCoefficientSize()),
        _nParam(evaluation.getEvaluator()->getParameterSize()),
        _nUnified(_nParam + _nCoeff),
        _penalty(0.0), 
        _parameters(lsst::ndarray::allocate(_nUnified)),
        _evaluation(evaluation)
    {
        useHybrid();
        setTol(fTol, gTol);
        setMinStep(minStep);
        setMaxIter(maxIter);
        setTau(tau);
    };

    typedef Eigen::Block<Eigen::VectorXd> VectorBlock;
    typedef Eigen::Block<Eigen::MatrixXd> MatrixBlock;

    virtual void calculateF(
        Eigen::VectorXd const & x, 
        Eigen::VectorXd & f
    ) const {
        lsst::ndarray::viewAsEigen(_parameters) = x;
        _penalty = _evaluation.getEvaluator()->clipToBounds(_parameters);

        _evaluation.update(
            _parameters[lsst::ndarray::view(0, _nParam)], 
            _parameters[lsst::ndarray::view(_nParam, _nUnified)]
        );

        f << lsst::ndarray::viewAsEigen(_evaluation.getResiduals());

        if (_penalty > 0.0) {
            f *= (1.0 + _penalty);
        }
    }

    virtual void calculateJ(
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & f,
        Eigen::MatrixXd & j
    ) const {
        j << lsst::ndarray::viewAsEigen(_evaluation.getResidualsJacobian()),
             lsst::ndarray::viewAsEigen(_evaluation.getModelMatrix());
        if (_penalty > 0.0) {
            j *= 1.0 + _penalty;
        }
    }

private:
    int _nCoeff;
    int _nParam;
    int _nUnified;
    mutable double _penalty;
    lsst::ndarray::Array<double,1,1> _parameters;
    lsst::meas::multifit::Evaluation & _evaluation;
};

}
namespace lsst {
namespace meas {
namespace multifit {


GaussianDistribution GaussNewtonOptimizer::solve(
    BaseEvaluator::Ptr const & evaluator,
    double const fTol, double const gTol, 
    double const minStep, 
    int const maxIter, 
    double const tau, 
    bool retryWithSvd
) {
    int nCoeff = evaluator->getCoefficientSize();
    int nParam = evaluator->getParameterSize();
    Eigen::VectorXd unified(nCoeff + nParam);
    Eigen::VectorXd residual(evaluator->getDataSize());
 
    Evaluation evaluation(evaluator);
    unified << ndarray::viewAsEigen(evaluation.getParameters()),
               ndarray::viewAsEigen(evaluation.getCoefficients());

    ::Solver solver(evaluation, fTol, gTol, minStep, maxIter, tau);     
    _solverSuccess = solver.solve(unified, residual);

    if(!_solverSuccess && retryWithSvd) {
        solver.useSVD();
        _solverSuccess = solver.solve(unified, residual);
    }
   
    Eigen::MatrixXd covariance;
    solver.getCovariance(covariance);
    return GaussianDistribution(unified, covariance);
}

}}}
