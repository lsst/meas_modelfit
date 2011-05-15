#include "lsst/meas/multifit/BruteForceSourceOptimizer.h"
#include "lsst/meas/multifit/Evaluation.h"

#include <Eigen/Cholesky>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

void BruteForceSourceOptimizer::solve(Evaluator::Ptr const & evaluator, int nTestPoints) {
    static double const SQRT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    int nParam = evaluator->getParameterSize();
    int nCoeff = evaluator->getCoefficientSize();
    if (_parameters.getSize<0>() != nTestPoints || _parameters.getSize<1>() != nParam) {
        _parameters = ndarray::allocate(nTestPoints, nParam);
    }
    if (_objectiveValues.getSize<0>() != nTestPoints) {
        _objectiveValues = ndarray::allocate(nTestPoints);
    }
    if (_bestCoefficients.getSize<0>() != nCoeff) {
        _bestCoefficients = ndarray::allocate(nCoeff);
        _coefficientCovariance = ndarray::allocate(nCoeff, nCoeff);
    }
    evaluator->writeInitialParameters(_parameters[0]);
    for (int i = 1; i < nTestPoints; ++i) {
        _parameters[i] = _parameters[0];
    }
    grid::Source const & source = evaluator->getGrid()->sources[0];
    afw::geom::ellipses::Ellipse psfEllipse = source.getLocalPsf()->computeMoments();
    afw::geom::ellipses::Quadrupole fullMoments(
        multifit::EllipseCore(
            source.object.getEllipticity()->getValue(),
            source.object.getRadius()->getValue()
        )
    );
    afw::geom::ellipses::Quadrupole psfMoments(psfEllipse.getCore());
    double minIxx = std::max(fullMoments.getIXX() - psfMoments.getIXX(), SQRT_EPS);
    double minIyy = std::max(fullMoments.getIYY() - psfMoments.getIYY(), SQRT_EPS);
    double minIxy = fullMoments.getIXY() - psfMoments.getIXY();
    minIxy = ((minIxy > 0) ? 1 : -1) * (SQRT_EPS + std::min(std::abs(minIxy), std::sqrt(minIxx * minIyy)));
    double dIxx = (fullMoments.getIXX() - minIxx) / nTestPoints;
    double dIyy = (fullMoments.getIYY() - minIyy) / nTestPoints;
    double dIxy = (fullMoments.getIXY() - minIxy) / nTestPoints;
    Evaluation evaluation(evaluator);
    _bestIndex = -1;
    double bestValue = std::numeric_limits<double>::infinity();
    for (int i = 0; i < nTestPoints; ++i) {
        afw::geom::ellipses::Quadrupole q(minIxx + i * dIxx, minIyy + i * dIyy, minIxy + i * dIxy);
        Ellipse ellipse(q, afw::geom::Point2D(source.object.getPosition()->getValue()));
        source.object.readEllipse(_parameters[i].getData(), ellipse);
        evaluation.update(_parameters[i]);
        _objectiveValues[i] = evaluation.getObjectiveValue();
        if (evaluation.getObjectiveValue() < bestValue) {
            _bestIndex = i;
            bestValue = evaluation.getObjectiveValue();
            _bestCoefficients.deep() = evaluation.getCoefficients();
            _coefficientCovariance.deep() = evaluation.getCoefficientFisherMatrix();
        }
    }
    Eigen::LDLT<Eigen::MatrixXd> ldlt(ndarray::viewAsTransposedEigen(_coefficientCovariance));
    ndarray::viewAsTransposedEigen(_coefficientCovariance).setIdentity();
    ldlt.solveInPlace(ndarray::viewAsTransposedEigen(_coefficientCovariance).setIdentity());
}


}}}
