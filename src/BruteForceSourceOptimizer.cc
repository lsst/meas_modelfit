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
    grid::Source const & ellipseSource = evaluator->getGrid()->sources[0];
    grid::Source const & psfSource = evaluator->getGrid()->sources[1];
    afw::geom::ellipses::Ellipse psfEllipse = ellipseSource.getLocalPsf()->computeMoments();
    afw::geom::ellipses::Quadrupole fullMoments(
        multifit::EllipseCore(
            ellipseSource.object.getEllipticity()->getValue(),
            ellipseSource.object.getRadius()->getValue()
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
    double dIxx = (fullMoments.getIXX() - minIxx) / (nTestPoints + 1);
    double dIyy = (fullMoments.getIYY() - minIyy) / (nTestPoints + 1);
    double dIxy = (fullMoments.getIXY() - minIxy) / (nTestPoints + 1);
    Evaluation evaluation(evaluator);
    _bestIndex = -1;
    double bestValue = std::numeric_limits<double>::infinity();
    _bestIsSafe = false;
    minIxx += dIxx;
    minIyy += dIyy;
    minIxy += dIxy;
    for (int i = 0; i < nTestPoints; ++i) {
        afw::geom::ellipses::Quadrupole q(minIxx + i * dIxx, minIyy + i * dIyy, minIxy + i * dIxy);
        Ellipse ellipse(q, afw::geom::Point2D(ellipseSource.object.getPosition()->getValue()));
        ellipseSource.object.readEllipse(_parameters[i].getData(), ellipse);
        try {
            evaluation.update(_parameters[i]);
            _objectiveValues[i] = evaluation.getObjectiveValue();
        } catch (...) {
            continue;
        }
        if (evaluation.getObjectiveValue() < bestValue) {
            double flux = ellipseSource.computeFluxMean(_parameters[i], evaluation.getCoefficients());
            flux += psfSource.computeFluxMean(_parameters[i], evaluation.getCoefficients());
            double condition = flux / ndarray::viewAsEigen(evaluation.getCoefficients()).norm();
            if (condition < 1E-10) {
                if (_bestIsSafe) continue;
            } else {
                _bestIsSafe = true;
            }
            _bestIndex = i;
            bestValue = evaluation.getObjectiveValue();
            _bestCoefficients.deep() = evaluation.getCoefficients();
            _coefficientCovariance.deep() = evaluation.getCoefficientFisherMatrix();
        }
    }
    if (_bestIndex >= 0) {
	Eigen::LDLT<Eigen::MatrixXd> ldlt(ndarray::viewAsTransposedEigen(_coefficientCovariance));
	ndarray::viewAsTransposedEigen(_coefficientCovariance).setIdentity();
	ldlt.solveInPlace(ndarray::viewAsTransposedEigen(_coefficientCovariance).setIdentity());
    }
    return _bestIndex >= 0;
}


}}}
