#include "lsst/meas/multifit/BruteForceSourceOptimizer.h"
#include "lsst/meas/multifit/Evaluation.h"

#include <Eigen/LU>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

GaussianDistribution::Ptr BruteForceSourceOptimizer::solve(Evaluator::Ptr const & evaluator, int n) {
    static double const SQRT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    if (_parameters.getSize<0>() != n || _parameters.getSize<1>() != evaluator->getParameterSize()) {
        _parameters = ndarray::allocate(n, evaluator->getParameterSize());
    }
    if (_objectiveValues.getSize<0>() != n) {
        _objectiveValues = ndarray::allocate(n);
    }
    evaluator->writeInitialParameters(_parameters[0]);
    for (int i = 1; i < n; ++i) {
        _parameters[i] = _parameters[0];
    }
    grid::Source const & source = evaluator->getGrid()->sources[0];
    afw::geom::ellipses::Ellipse psfEllipse = source.getLocalPsf()->computeMoments();
    afw::geom::ellipses::Quadrupole m(
        multifit::EllipseCore(
            source.object.getEllipticity()->getValue(),
            source.object.getRadius()->getValue()
        )
    );
    afw::geom::ellipses::Quadrupole p(psfEllipse.getCore());
    double minIxx = std::max(m.getIXX() - p.getIXX(), SQRT_EPS);
    double minIyy = std::max(m.getIYY() - p.getIYY(), SQRT_EPS);
    double minIxy = m.getIXY() - p.getIXY();
    minIxy = ((minIxy > 0) ? 1 : -1) * (SQRT_EPS + std::min(std::abs(minIxy), std::sqrt(minIxx * minIyy)));
    double dIxx = (m.getIXX() - minIxx) / n;
    double dIyy = (m.getIYY() - minIyy) / n;
    double dIxy = (m.getIXY() - minIxy) / n;
    Evaluation evaluation(evaluator);
    int bestIndex = 0;
    double bestValue = std::numeric_limits<double>::infinity();
    Eigen::VectorXd mu(evaluator->getParameterSize() + evaluator->getCoefficientSize());
    Eigen::MatrixXd fisherMatrix;
    for (int i = 0; i < n; ++i) {
        afw::geom::ellipses::Quadrupole q(minIxx + i * dIxx, minIyy + i * dIyy, minIxy + i * dIxy);
        Ellipse ellipse(q, afw::geom::Point2D(source.object.getPosition()->getValue()));
        source.object.readEllipse(_parameters[i].getData(), ellipse);
        evaluation.update(_parameters[i]);
        _objectiveValues[i] = evaluation.getObjectiveValue();
        if (evaluation.getObjectiveValue() < bestValue) {
            bestIndex = i;
            bestValue = evaluation.getObjectiveValue();
            mu.segment(0, evaluator->getParameterSize()) = ndarray::viewAsEigen(_parameters[i]);
            mu.segment(evaluator->getParameterSize(), evaluator->getCoefficientSize()) 
                = ndarray::viewAsEigen(evaluation.getCoefficients());
            fisherMatrix = ndarray::viewAsEigen(evaluation.getCoefficientFisherMatrix());
        }
    }
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(
        evaluator->getParameterSize() + evaluator->getCoefficientSize(),
        evaluator->getParameterSize() + evaluator->getCoefficientSize()
    );
    sigma.block(
        evaluator->getParameterSize(), evaluator->getParameterSize(), 
        evaluator->getCoefficientSize(), evaluator->getCoefficientSize()
    ) = sigma.inverse();
    return GaussianDistribution::Ptr(new GaussianDistribution(mu, sigma));
}


}}}
