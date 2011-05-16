#include "lsst/meas/multifit/BruteForceSourceOptimizer.h"
#include "lsst/meas/multifit/Evaluation.h"

#include <Eigen/Cholesky>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

bool BruteForceSourceOptimizer::solve(Evaluator::Ptr const & evaluator, int nTestPoints) {
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
    double minIxx = fullMoments.getIXX() - psfMoments.getIXX();
    double minIyy = fullMoments.getIYY() - psfMoments.getIYY();
    double minIxy = fullMoments.getIXY() - psfMoments.getIXY();
    if ((minIxx < SQRT_EPS) || (minIyy < SQRT_EPS) || (minIxy * minIxy > minIxx * minIyy)) {
        minIxx = SQRT_EPS;
        minIyy = SQRT_EPS;
        minIxy = 0.0;
    }
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
        try {
            evaluation.update(_parameters[i]);
            _objectiveValues[i] = evaluation.getObjectiveValue();
        } catch (...) {
            continue;
        }
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
    return _bestIndex >= 0;
}


}}}
