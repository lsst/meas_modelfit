#include "lsst/meas/multifit/GaussNewtonOptimizer.h"
#include "lsst/meas/algorithms/shapelet/NLSolver.h"

#include <Eigen/Core>
#include <cfloat>
#include "lsst/ndarray/eigen.h"


namespace lsst {
namespace meas {
namespace multifit {
namespace {

class Solver : public algorithms::shapelet::NLSolver {
public:
    Solver(Evaluation & evaluation, 
           double const fTol, double const gTol, 
           double const minStep, int const maxIter, 
           double const tau
    ) : _nCoeff(evalutaion.getEvaluator()->getCoefficientSize()),
        _nParam(evalutaion.getEvaluator()->getParameterSize()),
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
        ndarray::Array<double const 1, 1> u = ndarray::viewVectorAsArray(x);
        _evaluation.update(u[ndarray::view(0, _nParam)], u[ndarray::view(_nParam, x.size())]);

        f << ndarray::viewAsEigen(_evalutaion.getResiduals());
    }

    virtual void calculateJ(
        Eigen::VectorXd const & x,
        Eigen::VectorXd const & f,
        Eigen::MatrixXd & j
    ) const {
        j << ndarray::viewAsEigen(_evaluation.getResidualsJacobian()),
             ndarray::viewAsEigen(_evaluation.getModelMatrix());
    }

private:
    int _nCoeff, _nParam;
    Evaluation & _evaluation;
};

}

GuassianDistribution GaussNewtonOptimizer::solve(
    BaseEvaluator::Ptr const & evaluator,
    double const fTol, double const gTol, 
    double const minStep, 
    int const maxIter, 
    double const tau, 
    bool retryWithSvd
)
{
    int nCoeff = evaluator->getCoefficientSize();
    int nParam = evaluator->getParameterSize();
    Eigen::VectorXd unified(nCoeff + nParam);
    Eigen::VectorXd residual(evaluator->getDataSize());
 
    unified << ndarray::viewAsEigen(_evaluation.getParameters()),
               ndarray::viewAsEigen(_evaluation.getCoefficients());

    ::Solver solver(evaluation, fTol, gTol, minStep, maxIter, tau);     
    _solverSuccess = solver.solve(unified, residual);

    if(!_solverSuccess && retryWithSvd) {
        solver.useSVD();
        _solverSuccess = solver.solve(unified, residual);
    }
    
    return GaussianDistribution(unified, solver.getCovariance());
}

}}}
