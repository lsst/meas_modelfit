#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/meas/multifit/qp.h"
#include "lsst/ndarray/eigen.h"
#include <Eigen/Array>
#include <Eigen/SVD>

namespace lsst { namespace meas { namespace multifit {

namespace {

enum ProductEnum {
    MODEL_MATRIX = 0,
    MODEL_MATRIX_DERIVATIVE,
    COEFFICIENTS,
    RESIDUALS,
    RESIDUALS_JACOBIAN,
    MODEL_VECTOR,
    OBJECTIVE_VALUE,
    COEFFICIENT_FISHER_MATRIX,
    COEFFICIENT_COVARIANCE_MATRIX,
    FACTORIZATION,
    PRODUCT_COUNT
};

template <ProductEnum product>
struct Bit {
    static int const flag = 1 << product;
    static bool test(int status) { return flag & status; }
    static void set(int & status) { status |= flag; }
    static void reset(int & status) { status &= ~flag; }
};

static int const coefficient_dependencies = 
    Bit<MODEL_VECTOR>::flag |
    Bit<RESIDUALS>::flag |
    Bit<RESIDUALS_JACOBIAN>::flag |
    Bit<OBJECTIVE_VALUE>::flag
    ;

} // anonymous

class Evaluation::Factorization {
public:
    
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign | Eigen::RowMajor> MatrixRM;

    Eigen::JacobiSVD<MatrixRM> svd;
    int n1;
    int n2;
    int n;

    Eigen::VectorXd workspace;

    void factor(ndarray::Array<Pixel,2,1> const & modelMatrix, double svThreshold) {
        if (modelMatrix.getSize<0>() < modelMatrix.getSize<1>()) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::RuntimeErrorException,
                "Not enough data points to fit model."
            );
        }
        svd.compute(modelMatrix.asEigen(), Eigen::ComputeThinU | Eigen::ComputeThinV);
        n = n1 = svd.singularValues().size();
        workspace.resize(n);
        n2 = 0;
        svThreshold *= svd.singularValues()[0];
        while (n1 >= 1 && svd.singularValues()[n1 - 1] < svThreshold) {
            --n1;
            ++n2;
        }
    }

    void solve(
        ndarray::Array<Pixel,1,1> const & coefficients,
        ndarray::Array<Pixel const,1,1> const & data,
        ndarray::Array<Pixel const,2,2> const & constraintMatrix,
        ndarray::Array<Pixel const,1,1> const & constraintVector
    ) {
        Eigen::MatrixXd A = constraintMatrix.asEigen() * svd.matrixV().block(0, 0, n, n1);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n1, n1);
        G.diagonal() = svd.singularValues().segment(0, n1);
        Eigen::VectorXd c = - svd.matrixU().block(0, 0, svd.matrixU().rows(), n1).transpose() 
            * data.asEigen();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n1);
        QPSolver(G, c).inequality(A, constraintVector.asEigen()).solve(x);
        coefficients.asEigen() = svd.matrixV().block(0, 0, n, n1) * x;
    }

    void fillFisherMatrix(ndarray::Array<Pixel,2,2> const & matrix) {
        Eigen::MatrixXd tmp = svd.singularValues().segment(0, n1).asDiagonal() 
            * svd.matrixV().block(0, 0, n, n1).transpose();
        matrix.asEigen().selfAdjointView<Eigen::Upper>().rankUpdate(tmp.adjoint());
        matrix.asEigen().triangularView<Eigen::StrictlyLower>() =
            matrix.asEigen().triangularView<Eigen::StrictlyUpper>().transpose();
    }

    void fillCovarianceMatrix(ndarray::Array<Pixel,2,2> const & matrix) {
        Eigen::MatrixXd tmp 
            = svd.singularValues().segment(0, n1).array().inverse().matrix().asDiagonal() 
            * svd.matrixV().block(0, 0, n, n1).transpose();
        matrix.asEigen().selfAdjointView<Eigen::Upper>().rankUpdate(tmp.adjoint());
        matrix.asEigen().triangularView<Eigen::StrictlyLower>() =
            matrix.asEigen().triangularView<Eigen::StrictlyUpper>().transpose();
    }

};

Evaluation::Evaluation(BaseEvaluator::Ptr const & evaluator, double const svThreshold) : 
    _products(0), _evaluator(evaluator), _factorization(new Factorization),
    _parameters(ndarray::allocate(_evaluator->getParameterCount())),
    _svThreshold(svThreshold)
{
    _evaluator->writeInitialParameters(_parameters);
    initialize();
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator,
    lsst::ndarray::Array<Pixel const,1,1> const & parameters, 
    double const svThreshold
) : 
    _products(0), _evaluator(evaluator), _factorization(new Factorization),
    _parameters(ndarray::copy(parameters)),
    _svThreshold(svThreshold)
{
    initialize();
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator,
    Eigen::VectorXd const & parameters, 
    double const svThreshold
) : 
    _products(0), _evaluator(evaluator),  _factorization(new Factorization),
    _parameters(ndarray::allocate(parameters.size())),
    _svThreshold(svThreshold)
{
    _parameters.asEigen() = parameters;
    initialize();
}

void Evaluation::update(lsst::ndarray::Array<double const,1,1> const & parameters) {
    assert(parameters.getSize<0>() == _parameters.getSize<0>());
    _parameters.deep() = parameters;
    _products = 0;
}

void Evaluation::update(Eigen::VectorXd const & parameters) {
    assert(parameters.size() == _parameters.getSize<0>());
    _parameters.asEigen() = parameters;
    _products = 0;
}
    
void Evaluation::update(
    lsst::ndarray::Array<double const,1,1> const & parameters, 
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients
) {
    update(parameters);
    setCoefficients(coefficients);
}

void Evaluation::update(
    Eigen::VectorXd const & parameters,
    Eigen::VectorXd const & coefficients
) {
    update(parameters);
    setCoefficients(coefficients);
}

void Evaluation::setCoefficients(lsst::ndarray::Array<Pixel const,1,1> const & coefficients) {
    assert(coefficients.size() == _evaluator->getCoefficientCount());
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientCount());
    }
    _coefficients.deep() = coefficients;
    _products &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::set(_products);
}

void Evaluation::setCoefficients(Eigen::VectorXd const & coefficients) {
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientCount());
    }
    assert(coefficients.size() == _evaluator->getCoefficientCount());
    _coefficients.asEigen() = coefficients;
    _products &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::set(_products);
}

void Evaluation::solveCoefficients() {
    _products &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::reset(_products);
    ensureCoefficients();
}

void Evaluation::ensureModelMatrix() const {
    if (Bit<MODEL_MATRIX>::test(_products)) return;
    if (_modelMatrix.getData() == 0) {
        _modelMatrix = ndarray::allocate(_evaluator->getPixelCount(), _evaluator->getCoefficientCount());
    }
    _evaluator->_evaluateModelMatrix(_modelMatrix, _parameters);
    Bit<MODEL_MATRIX>::set(_products);
}

void Evaluation::ensureModelMatrixDerivative() const {
    if (Bit<MODEL_MATRIX_DERIVATIVE>::test(_products)) return;
    ensureModelMatrix();
    if (_modelMatrixDerivative.getData() == 0) {
        _modelMatrixDerivative = ndarray::allocate(
            _evaluator->getParameterCount(), _evaluator->getPixelCount(), _evaluator->getCoefficientCount()
        );
    }
    _evaluator->_evaluateModelMatrixDerivative(_modelMatrixDerivative, _modelMatrix, _parameters);
    Bit<MODEL_MATRIX_DERIVATIVE>::set(_products);
}

void Evaluation::ensureCoefficients() const {
    if (Bit<COEFFICIENTS>::test(_products)) return;
    ensureFactorization();
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientCount());
    }
    _factorization->solve(
        _coefficients, _evaluator->getDataVector(),
        _evaluator->getConstraintMatrix(), _evaluator->getConstraintVector()
    );
    Bit<COEFFICIENTS>::set(_products);
}

void Evaluation::ensureModelVector() const {
    if (Bit<MODEL_VECTOR>::test(_products)) return;
    ensureCoefficients();
    ensureModelMatrix();
    if (_modelVector.getData() == 0) {
        _modelVector = ndarray::allocate(_evaluator->getPixelCount());
    }
    _modelVector.asEigen() 
        = _modelMatrix.asEigen() * _coefficients.asEigen();
    Bit<MODEL_VECTOR>::set(_products);
}

void Evaluation::ensureResiduals() const {
    if (Bit<RESIDUALS>::test(_products)) return;
    ensureModelVector();
    if (_residuals.getData() == 0) {
        _residuals = ndarray::allocate(_evaluator->getPixelCount());
    }
    _residuals.asEigen() = _modelVector.asEigen() 
        - _evaluator->getDataVector().asEigen();
    Bit<RESIDUALS>::set(_products);
}

void Evaluation::ensureResidualsJacobian() const {
    if (Bit<RESIDUALS_JACOBIAN>::test(_products)) return;
    ensureCoefficients();
    ensureModelMatrixDerivative();
    if (_residualsJacobian.getData() == 0) {
        _residualsJacobian = ndarray::allocate(
            _evaluator->getPixelCount(), _evaluator->getParameterCount()
        );
    }
    ndarray::EigenView<Pixel,1,1> coefficients(_coefficients);
    ndarray::EigenView<Pixel,2,2> jac(_residualsJacobian);
    for (int n = 0; n < jac.cols(); ++n) {
        jac.col(n) = _modelMatrixDerivative[n].asEigen() * coefficients;
    }
    Bit<RESIDUALS_JACOBIAN>::set(_products);
}

void Evaluation::ensureCoefficientFisherMatrix() const {
    if (Bit<COEFFICIENT_FISHER_MATRIX>::test(_products)) return;
    if (_coefficientFisherMatrix.getData() == 0) {
        _coefficientFisherMatrix = ndarray::allocate(
            _evaluator->getCoefficientCount(), _evaluator->getCoefficientCount()
        );
    }
    if (Bit<FACTORIZATION>::test(_products)) {
        _factorization->fillFisherMatrix(_coefficientFisherMatrix);
    } else {
        ensureModelMatrix();
        _coefficientFisherMatrix.asEigen().part<Eigen::SelfAdjoint>()
            = _modelMatrix.asEigen().transpose()
            * _modelMatrix.asEigen();
    }
    Bit<COEFFICIENT_FISHER_MATRIX>::set(_products);
}

void Evaluation::ensureCoefficientCovarianceMatrix() const {
    if (Bit<COEFFICIENT_COVARIANCE_MATRIX>::test(_products)) return;
    ensureFactorization();
    if (_coefficientCovarianceMatrix.getData() == 0) {
        _coefficientCovarianceMatrix = ndarray::allocate(
            _evaluator->getCoefficientCount(), _evaluator->getCoefficientCount()
        );
    }
    _factorization->fillCovarianceMatrix(_coefficientCovarianceMatrix);
    Bit<COEFFICIENT_COVARIANCE_MATRIX>::set(_products);
}

void Evaluation::ensureObjectiveValue() const {
    if (Bit<OBJECTIVE_VALUE>::test(_products)) return;
    ensureResiduals();
    _objectiveValue = 0.5 * _residuals.asEigen().squaredNorm();
    Bit<OBJECTIVE_VALUE>::set(_products);
}

void Evaluation::ensureFactorization() const {
    if (Bit<FACTORIZATION>::test(_products)) return;
    ensureModelMatrix();
    _factorization->factor(_modelMatrix, _svThreshold);
    Bit<FACTORIZATION>::set(_products);
}

void Evaluation::initialize() {
    if (_evaluator->getParameterCount() != _parameters.getSize<0>()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            (boost::format("Evaluator parameter size (%d) does not match parameter vector size (%d).")
             % _evaluator->getParameterCount() % _parameters.getSize<0>()).str()
        );
    }
}

Evaluation::~Evaluation() {}


}}} // namespace lsst::meas::multifit
