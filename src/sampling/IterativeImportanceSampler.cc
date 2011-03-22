#include "lsst/meas/multifit/sampling/IterativeImportanceSampler.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

#include <Eigen/Cholesky>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

double IterativeImportanceSampler::computeNormalizedPerplexity() const {
    Eigen::VectorXd w = ndarray::viewAsEigen(_samples.back().weight) / _samples.back().getWeightSum();
    return -(w.cwise() * w.cwise().log()).sum();
}

double IterativeImportanceSampler::computeEffectiveSampleSize() const {
    return 1.0 / (ndarray::viewAsEigen(_samples.back().weight) 
                  / _samples.back().getWeightSum()).cwise().square().sum();
}

void IterativeImportanceSampler::run(int size) {
    ndarray::Array<Pixel,2,2> modelArray =
        ndarray::allocate(_evaluator->getDataSize(), _evaluator->getCoefficientSize());
    _samples.push_back(
        Table::allocate(size, _evaluator->getParameterSize(), FISHER_LLT, _evaluator->getCoefficientSize())
    );
    ndarray::EigenView<Pixel,2,2> modelMatrix(modelArray);
    Eigen::MatrixXd fisher(_evaluator->getCoefficientSize(), _evaluator->getCoefficientSize());
    Eigen::LLT<Eigen::MatrixXd> llt;
    Table const & table = _samples.back();
    _importance.draw(table, _randomEngine);
    double sumDataLn2Pi = 0.5 * _evaluator->getDataSize() * std::log(2.0 * M_PI);
    double sumCoeffLn2Pi = 0.5 * _evaluator->getCoefficientSize() * std::log(2.0 * M_PI);
    for (int n = 0; n < size; ++n) {
        Record record = table[n];
        _evaluator->evaluateModelMatrix(modelArray, record.parameters);
        fisher.part<Eigen::SelfAdjoint>() = modelMatrix.transpose() * modelMatrix;
        record.nested.matrix.part<Eigen::LowerTriangular>() = llt.matrixL();
        record.nested.vector = 
            (modelMatrix.transpose() * ndarray::viewAsEigen(_evaluator->getDataVector())).lazy();
        llt.solveInPlace(record.nested.vector);
        record.nested.scalar = 0.5 * 
            (ndarray::viewAsEigen(_evaluator->getDataVector())
             - modelMatrix * record.nested.vector).squaredNorm();
        record.nested.scalar += 0.5 * _evaluator->getLogVarianceSum() + sumDataLn2Pi;
        /// target = scalar + 0.5 * ln |F/2pi| (likelihood marginalized with flat prior)
        record.target = record.nested.scalar 
            + record.nested.matrix.diagonal().cwise().log().sum() - sumCoeffLn2Pi;
    }
    table.computeWeights();
    _importance.update(table);
}

IterativeImportanceSampler::IterativeImportanceSampler(
    BaseEvaluator::Ptr const & evaluator,
    MixtureDistribution const & importance,
    SampleList const & samples,
    RandomEngine const & randomEngine
) :
    _evaluator(evaluator),
    _randomEngine(randomEngine),
    _importance(importance),
    _samples(samples)
{}

}}}} // namespace lsst::meas::multifit::sampling
