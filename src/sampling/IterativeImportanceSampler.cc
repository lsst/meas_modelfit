#include "lsst/meas/multifit/sampling/IterativeImportanceSampler.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

#include <Eigen/SVD>
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
    ndarray::EigenView<Pixel const,1,1> dataVector(_evaluator->getDataVector());
#ifdef USE_CHOLESKY
    Eigen::MatrixXd fisher(_evaluator->getCoefficientSize(), _evaluator->getCoefficientSize());
    Eigen::LLT<Eigen::MatrixXd> llt;
#else
    Eigen::SVD<Eigen::MatrixXd> svd;
#endif
    Table const & table = _samples.back();
    _importance.draw(table, _randomEngine);
    //double sumDataLn2Pi = 0.5 * _evaluator->getDataSize() * std::log(2.0 * M_PI);
    //double sumCoeffLn2Pi = 0.5 * _evaluator->getCoefficientSize() * std::log(2.0 * M_PI);
    for (int n = 0; n < size; ++n) {
        Record record = table[n];
        _evaluator->evaluateModelMatrix(modelArray, record.parameters);
#ifdef USE_CHOLESKY
        fisher.part<Eigen::SelfAdjoint>() = modelMatrix.transpose() * modelMatrix;
        llt.compute(fisher);
        record.nested.matrix.part<Eigen::LowerTriangular>() = llt.matrixL();
        record.nested.vector = (modelMatrix.transpose() * dataVector).lazy();
        llt.solveInPlace(record.nested.vector);
        double q = record.nested.matrix.diagonal().cwise().log().sum();
#else
        svd.compute(modelMatrix);
        record.nested.matrix.part<Eigen::SelfAdjoint>() = modelMatrix.transpose() * modelMatrix;
        svd.sort();
        Eigen::VectorXd sv = svd.singularValues();
        double q = 0.5 * sv.cwise().log().sum();
        for (int n = 0; n < sv.size(); ++n) {
            sv[n] = (sv[n] < sv[0] * 1E-8) ? 0.0 : 1.0 / sv[n];
        }
        Eigen::VectorXd z = svd.matrixU().transpose() * dataVector;
        record.nested.vector = svd.matrixV() * (sv.cwise() * z);
#endif
        record.nested.scalar = 0.5 * (dataVector - modelMatrix * record.nested.vector).squaredNorm();
        // record.nested.scalar += 0.5 * _evaluator->getLogVarianceSum() + sumDataLn2Pi;
        /// target = scalar + 0.5 * ln |F/2pi| (likelihood marginalized with flat prior)
        record.target = record.nested.scalar + q; // - sumCoeffLn2Pi;
    }
    table.computeWeights();
    _importance.update(table);
}

IterativeImportanceSampler::IterativeImportanceSampler(
    BaseEvaluator::Ptr const & evaluator,
    MixtureDistribution const & importance,
    RandomEngine const & randomEngine
) :
    _evaluator(evaluator),
    _randomEngine(randomEngine),
    _importance(importance)
{}

}}}} // namespace lsst::meas::multifit::sampling
