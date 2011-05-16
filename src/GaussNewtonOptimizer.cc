#include "lsst/meas/multifit/GaussNewtonOptimizer.h"
#include "lsst/meas/algorithms/shapelet/NLSolver.h"
#include "lsst/meas/multifit/Evaluation.h"

#include <Eigen/Core>
#include <cfloat>
#include "lsst/ndarray/eigen.h"

#include <iostream>

namespace {

class Solver : public lsst::meas::algorithms::shapelet::NLSolver {
public:
    Solver(lsst::meas::multifit::Evaluation & evaluation, 
           std::list<lsst::ndarray::Array<double, 1, 1> > & parameterPoints,
           double const fTol, double const gTol, 
           double const minStep, int const maxIter, 
           double const tau
    ) : _nCoeff(evaluation.getEvaluator()->getCoefficientSize()),
        _nParam(evaluation.getEvaluator()->getParameterSize()),
        _nUnified(_nParam + _nCoeff),
        _penalty(0.0), 
        _parameterPoints(parameterPoints),
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
        lsst::ndarray::Array<double, 1, 1> parameters = lsst::ndarray::copy(
            lsst::ndarray::viewVectorAsArray(x)
        );
        _parameterPoints.push_back(parameters);
        _penalty = _evaluation.getEvaluator()->clipToBounds(parameters);

        _evaluation.update(
            parameters[lsst::ndarray::view(0, _nParam)], 
            parameters[lsst::ndarray::view(_nParam, _nUnified)]
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
        //std::cerr << "In calculateJ!\n";
        j << lsst::ndarray::viewAsEigen(_evaluation.getResidualsJacobian()),
             lsst::ndarray::viewAsEigen(_evaluation.getModelMatrix());
        if (_penalty > 0.0) {
            j *= 1.0 + _penalty;
        }
        //std::cerr << j << "\n";
    }

private:
    int _nCoeff;
    int _nParam;
    int _nUnified;
    mutable double _penalty;
    std::list<lsst::ndarray::Array<double,1,1> > &_parameterPoints;
    lsst::meas::multifit::Evaluation & _evaluation;
};

}
namespace lsst {
namespace meas {
namespace multifit {

GaussianDistribution::Ptr GaussNewtonOptimizer::solve(
    BaseEvaluator::Ptr const & evaluator,
    double const fTol, double const gTol, 
    double const minStep, 
    int const maxIter, 
    double const tau, 
    bool const retryWithSvd
) {
    int nCoeff = evaluator->getCoefficientSize();
    int nParam = evaluator->getParameterSize();
    
    if (nCoeff + nParam > evaluator->getDataSize()) {
        _didConverge=false;
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Have fewer pixels than parameters. System is underdetermined"
        );
    }
    _parameterPoints.clear();
    Eigen::VectorXd unified(nCoeff + nParam);
    Eigen::MatrixXd covariance(nParam+nCoeff, nParam+nCoeff);
    Evaluation evaluation(evaluator);
    if (nParam == 0) {
        _didConverge = true;
        unified << ndarray::viewAsEigen(evaluation.getCoefficients());
        _parameterPoints.push_back(ndarray::copy(evaluation.getCoefficients()));
        covariance << ndarray::viewAsEigen(evaluation.getCoefficientFisherMatrix()).inverse();
        return GaussianDistribution::Ptr(
            new GaussianDistribution(unified, covariance)
        );
    }

    Eigen::VectorXd residual(evaluator->getDataSize());
 

    unified << ndarray::viewAsEigen(evaluation.getParameters()),
                ndarray::viewAsEigen(evaluation.getCoefficients());
    
    _parameterPoints.push_back(ndarray::copy(ndarray::viewVectorAsArray(unified)));

    ::Solver solver(evaluation, _parameterPoints, fTol, gTol, minStep, maxIter, tau);     
    _didConverge = solver.solve(unified, residual);

    if(!_didConverge && retryWithSvd) {
        solver.useSVD();
        _didConverge = solver.solve(unified, residual);
    }
   
    solver.getCovariance(covariance);
    return GaussianDistribution::Ptr(
        new GaussianDistribution(unified, covariance)
    );
}

lsst::ndarray::Array<const double, 2, 2> GaussNewtonOptimizer::getParameterPoints() const {
    if (_parameterPoints.empty()) {
        return lsst::ndarray::Array<const double, 2, 2>();
    }
    std::list<ndarray::Array<double, 1, 1> >::const_iterator i(_parameterPoints.begin());
    ndarray::Array<double, 2, 2> parameterPoints = ndarray::allocate(
        static_cast<int>(_parameterPoints.size()),
        i->getSize<0>()
    );

    for(int j=0 ; i != _parameterPoints.end(); ++j) {
        parameterPoints[j].deep() = *i;
    }
    return parameterPoints;
}

}}}
