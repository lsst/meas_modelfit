#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/meas/multifit/qp.h"
#include "lsst/ndarray/eigen.h"
#include <Eigen/Cholesky>
#include <Eigen/Array>
#include <Eigen/QR>
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

class Evaluation::LinearSolver : private boost::noncopyable {
public:

    /**
     *  Solve a linear-least squares-like problem for the coefficients x that minimize:
     *  @f[
     *      (y - A x)^T (y - A x)
     *  @f]
     *  where @f$y@f$ is the data vector and @f$A@f$ is the model matrix.
     *
     *  @param[in]     modelMatrix   Matrix @f$A@f$.
     *  @param[in]     fisherMatrix  Matrix @f$F = A^T A@f$.
     *  @param[in]     data          Vector @f$y@f$.
     *  @param[out]    coefficients  Vector @f$x@f$.
     */
    virtual void solve(
        ndarray::EigenView<Pixel const,2,2> const & modelMatrix,
        ndarray::EigenView<Pixel const,2,2> const & fisherMatrix,
        ndarray::EigenView<Pixel const,1,1> const & data,
        ndarray::EigenView<Pixel,1,1> coefficients
    ) = 0;

    virtual ~LinearSolver() {}
};

class Evaluation::CholeskySolver : public Evaluation::LinearSolver {
public:

    virtual void solve(
        ndarray::EigenView<Pixel const,2,2> const & modelMatrix,
        ndarray::EigenView<Pixel const,2,2> const & fisherMatrix,
        ndarray::EigenView<Pixel const,1,1> const & data,
        ndarray::EigenView<Pixel,1,1> coefficients
    ) {
        coefficients = (modelMatrix.transpose() * data).lazy();
        Eigen::LDLT<MatrixRM> ldlt(fisherMatrix);
        ldlt.solveInPlace(coefficients);
    }

};

class Evaluation::ConstrainedSolver : public Evaluation::LinearSolver {
public:

    virtual void solve(
        ndarray::EigenView<Pixel const,2,2> const & modelMatrix,
        ndarray::EigenView<Pixel const,2,2> const & fisherMatrix,
        ndarray::EigenView<Pixel const,1,1> const & data,
        ndarray::EigenView<Pixel,1,1> coefficients
    ) {
        ndarray::viewAsEigen(_rhs) = -(modelMatrix.transpose() * data).lazy();
        QPSolver(fisherMatrix.getArray(), _rhs).inequality(_matrix, _vector).solve(coefficients.getArray());
    }

    ConstrainedSolver(
        ndarray::Array<Pixel const,2,1> const & matrix,
        ndarray::Array<Pixel const,1,1> const & vector
    ) : _rhs(ndarray::allocate(matrix.getSize<1>())), _matrix(matrix), _vector(vector) {}

private:
    ndarray::Array<Pixel,1,1> _rhs;
    ndarray::Array<Pixel const,2,1> _matrix;
    ndarray::Array<Pixel const,1,1> _vector;
};

Evaluation::Evaluation(BaseEvaluator::Ptr const & evaluator, double const svThreshold) : 
    _products(0), _evaluator(evaluator), 
    _parameters(ndarray::allocate(_evaluator->getParameterSize())),
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
    _products(0), _evaluator(evaluator), 
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
    _products(0), _evaluator(evaluator), 
    _parameters(ndarray::allocate(parameters.size())),
    _svThreshold(svThreshold)
{
    ndarray::viewAsEigen(_parameters) = parameters;
    initialize();
}

void Evaluation::update(lsst::ndarray::Array<double const,1,1> const & parameters) {
    assert(parameters.getSize<0>() == _parameters.getSize<0>());
    _parameters.deep() = parameters;
    _products = 0;
}
void Evaluation::update(Eigen::VectorXd const & parameters) {
    assert(parameters.size() == _parameters.getSize<0>());
    ndarray::viewAsEigen(_parameters) = parameters;
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
    assert(coefficients.size() == _evaluator->getCoefficientSize());
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientSize());
    }
    _coefficients.deep() = coefficients;
    _products &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::set(_products);
}

void Evaluation::setCoefficients(Eigen::VectorXd const & coefficients) {
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientSize());
    }
    assert(coefficients.size() == _evaluator->getCoefficientSize());
    ndarray::viewAsEigen(_coefficients) = coefficients;
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
        _modelMatrix = ndarray::allocate(_evaluator->getDataSize(), _evaluator->getCoefficientSize());
    }
    _evaluator->_evaluateModelMatrix(_modelMatrix, _parameters);
    Bit<MODEL_MATRIX>::set(_products);
}

void Evaluation::ensureModelMatrixDerivative() const {
    if (Bit<MODEL_MATRIX_DERIVATIVE>::test(_products)) return;
    ensureModelMatrix();
    if (_modelMatrixDerivative.getData() == 0) {
        _modelMatrixDerivative = ndarray::allocate(
            _evaluator->getParameterSize(), _evaluator->getDataSize(), _evaluator->getCoefficientSize()
        );
    }
    _evaluator->_evaluateModelMatrixDerivative(_modelMatrixDerivative, _modelMatrix, _parameters);
    Bit<MODEL_MATRIX_DERIVATIVE>::set(_products);
}

void Evaluation::ensureCoefficients() const {
    if (Bit<COEFFICIENTS>::test(_products)) return;
    ensureCoefficientFisherMatrix();
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientSize());
    }
    ndarray::EigenView<Pixel const,2,2> modelMatrix(_modelMatrix);
    ndarray::EigenView<Pixel const,2,2> fisherMatrix(_coefficientFisherMatrix);
    ndarray::EigenView<Pixel const,1,1> data(_evaluator->getDataVector());
    ndarray::EigenView<Pixel,1,1> coeff(_coefficients);
 
    _svd.reset();
    //try {
        _solver->solve(modelMatrix, fisherMatrix, data, coeff);
        /*
    } catch (...) {
        
        _svd.reset(new Eigen::SVD<MatrixRM>(modelMatrix.svd()));
        //Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> U = svd.matrixU();
        //Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> V = svd.matrixV();
        //Eigen::Matrix<Pixel, Eigen::Dynamic, 1> s = svd.singularValues();
        Eigen::Matrix<Pixel, Eigen::Dynamic, 1, Eigen::AutoAlign> tmp;

        tmp = _svd->matrixU().transpose()*(data);
        for(int i=0; i < tmp.size(); ++i) {
            double s = _svd->singularValues()[i];
            if(s > _svThreshold) {
                tmp[i] /= s;
            } else {
                tmp[i] = 0.;
            }
        }
        coeff = _svd->matrixV()*tmp;
    }
        */
    Bit<COEFFICIENTS>::set(_products);
}

void Evaluation::ensureModelVector() const {
    if (Bit<MODEL_VECTOR>::test(_products)) return;
    ensureCoefficients();
    ensureModelMatrix();
    if (_modelVector.getData() == 0) {
        _modelVector = ndarray::allocate(_evaluator->getDataSize());
    }
    ndarray::viewAsEigen(_modelVector) 
        = ndarray::viewAsEigen(_modelMatrix) * ndarray::viewAsEigen(_coefficients);
    Bit<MODEL_VECTOR>::set(_products);
}

void Evaluation::ensureResiduals() const {
    if (Bit<RESIDUALS>::test(_products)) return;
    ensureModelVector();
    if (_residuals.getData() == 0) {
        _residuals = ndarray::allocate(_evaluator->getDataSize());
    }
    ndarray::viewAsEigen(_residuals) = ndarray::viewAsEigen(_modelVector) 
        - ndarray::viewAsEigen(_evaluator->getDataVector());
    Bit<RESIDUALS>::set(_products);
}

void Evaluation::ensureResidualsJacobian() const {
    if (Bit<RESIDUALS_JACOBIAN>::test(_products)) return;
    ensureCoefficients();
    ensureModelMatrixDerivative();
    if (_residualsJacobian.getData() == 0) {
        _residualsJacobian = ndarray::allocate(
            _evaluator->getDataSize(), _evaluator->getParameterSize()
        );
    }
    ndarray::EigenView<Pixel,1,1> coeffVec(_coefficients);
    ndarray::EigenView<Pixel,2,2> jac(_residualsJacobian);
    for (int n = 0; n < jac.cols(); ++n) {
        jac.col(n) = ndarray::viewAsEigen(_modelMatrixDerivative[n]) * coeffVec;
    }
    Bit<RESIDUALS_JACOBIAN>::set(_products);
}

void Evaluation::ensureObjectiveValue() const {
    if (Bit<OBJECTIVE_VALUE>::test(_products)) return;
    ensureResiduals();
    _objectiveValue = 0.5 * ndarray::viewAsEigen(_residuals).squaredNorm();
    Bit<OBJECTIVE_VALUE>::set(_products);
}

void Evaluation::ensureCoefficientFisherMatrix() const {
    if (Bit<COEFFICIENT_FISHER_MATRIX>::test(_products)) return;
    ensureModelMatrix();
    if (_coefficientFisherMatrix.getData() == 0) {
        _coefficientFisherMatrix = ndarray::allocate(
            _evaluator->getCoefficientSize(), _evaluator->getCoefficientSize()
        );
    }
    ndarray::viewAsEigen(_coefficientFisherMatrix).part<Eigen::SelfAdjoint>()
        = ndarray::viewAsTransposedEigen(_modelMatrix)
        * ndarray::viewAsEigen(_modelMatrix);
    Bit<COEFFICIENT_FISHER_MATRIX>::set(_products);
}

void Evaluation::initialize() {
    if (_evaluator->getParameterSize() != _parameters.getSize<0>()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            (boost::format("Evaluator parameter size (%d) does not match parameter vector size (%d).")
             % _evaluator->getParameterSize() % _parameters.getSize<0>()).str()
        );
    }
    if (_evaluator->getConstraintSize() > 0) {
        _solver.reset(
            new ConstrainedSolver(_evaluator->getConstraintMatrix(), _evaluator->getConstraintVector())
        );
    } else {
        _solver.reset(new CholeskySolver());
    }
}

Evaluation::~Evaluation() {}


}}} // namespace lsst::meas::multifit
