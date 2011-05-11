#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/ndarray/eigen.h"
#include <Eigen/Cholesky>
#include <Eigen/Array>
#include <Eigen/QR>

namespace lsst { namespace meas { namespace multifit {

namespace {

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign | Eigen::RowMajor> MatrixRM;

enum ProductEnum {
    MODEL_MATRIX = 0,
    MODEL_MATRIX_DERIVATIVE,
    COEFFICIENTS,
    RESIDUALS,
    RESIDUALS_JACOBIAN,
    OBJECTIVE_VALUE,
    COEFFICIENT_FISHER_MATRIX,
    COEFFICIENT_FISHER_FACTOR,
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
     *      (y - A x)^T (y - A x) + (x - \mu)^T \Sigma^{-1} (x - \mu)^T
     *  @f]
     *  where @f$y@f$ is the data vector, @f$A@f$ is the model matrix, and 
     *  @f$\mu@f$ and @f$\Sigma@f$ are the mean and covariance
     *  of a Gaussian prior.
     *
     *  @param[in]     modelMatrix   Matrix @f$A@f$.  Always provided.
     *  @param[in]     data          Vector @f$y@f$.  Always provided.
     *  @param[out]    coefficients  Vector @f$x@f$.  Must always be filled.
     *  @param[in,out] fisherMatrix  Matrix @f$F = A^T A + \Sigma^{-1}@f$.  May or may not be provided,
     *                               and may or may not be filled.
     *  @param[in,out] fisherFactor  Lower triangular Cholesky factor of @f$F@f$.  May or may not be provided,
     *                               and may or may not be filled.
     *  @param[in,out] status        Bitflag indicating whether fisherMatrix and fisherFactor are valid
     *                               on input and output.
     */
    virtual void solve(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double const,1,1> const & data,
        ndarray::EigenView<double,1,1> coefficients,
        ndarray::Array<double,2,2> & fisherMatrix,
        ndarray::Array<double,2,2> & fisherFactor,
        int & status
    ) = 0;

    /**
     *  Compute the Lower-triangular Cholesky factor of the Fisher matrix @f$F@f$ as defined in solve().
     *
     *  @param[in]     modelMatrix   Matrix @f$A@f$.  Always provided.
     *  @param[out]    fisherFactor  Lower triangular Cholesky factor of @f$F@f$.  Must always be filled.
     *  @param[in,out] fisherMatrix  Matrix @f$F = A^T A + \Sigma^{-1}@f$.  May or may not be provided,
     *                               and may or may not be filled.
     *  @param[in,out] status        Bitflag indicating whether fisherMatrix is valid
     *                               on input and output.
     */
    virtual void computeFisherFactor(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double,2,2> fisherFactor,
        ndarray::Array<double,2,2> & fisherMatrix,
        int & status
    ) = 0;

    /**
     *  Compute the Fisher matrix @f$F@f$ as defined in solve().
     *
     *  @param[in]     modelMatrix   Matrix @f$A@f$.  Always provided.
     *  @param[out]    fisherMatrix  Matrix @f$F = A^T A + \Sigma^{-1}@f$.  Must always be filled.
     *  @param[in,out] fisherFactor  Lower triangular Cholesky factor of @f$F@f$.  May or may not be
     *                               provided, and may or may not be filled.
     *  @param[in,out] status        Bitflag indicating whether fisherFactor is valid
     *                               on input and output.
     */
    virtual void computeFisherMatrix(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double,2,2> fisherFactor,
        ndarray::Array<double,2,2> & fisherMatrix,
        int & status
    ) = 0;

    /**
     *  Set the parameters of the Gaussian prior.
     */
    virtual void setPrior(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma) {}

    /**
     *  Set the location parameter of the Gaussian prior, leaving the scale unchanged.
     */
    virtual void setPrior(Eigen::VectorXd const & mu) {}

    virtual ~LinearSolver() {}
};

class Evaluation::CholeskySolver : public Evaluation::LinearSolver {
public:
    typedef Eigen::LLT<MatrixRM> LLT;

    virtual void solve(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double const,1,1> const & data,
        ndarray::EigenView<double,1,1> coefficients,
        ndarray::Array<double,2,2> & fisherMatrix,
        ndarray::Array<double,2,2> & fisherFactor,
        int & status
    ) {
        if (!Bit<COEFFICIENT_FISHER_FACTOR>::test(status)) {
            if (fisherFactor.getData() == 0)
                fisherFactor = ndarray::allocate(modelMatrix.cols(), modelMatrix.cols());
            computeFisherFactor(modelMatrix, ndarray::viewAsEigen(fisherFactor), fisherMatrix, status);
            Bit<COEFFICIENT_FISHER_FACTOR>::set(status);
        }
        _rhs = (modelMatrix.transpose() * data).lazy();
        _rhs += _priorVector;
        ndarray::viewAsEigen(fisherFactor).part<Eigen::LowerTriangular>()
            .solveTriangularInPlace(_rhs);
        ndarray::viewAsEigen(fisherFactor).part<Eigen::LowerTriangular>().transpose()
            .solveTriangularInPlace(_rhs);
        coefficients = _rhs;
    }

    virtual void computeFisherFactor(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double,2,2> fisherFactor,
        ndarray::Array<double,2,2> & fisherMatrix,
        int & status
    ) {
        if (!Bit<COEFFICIENT_FISHER_MATRIX>::test(status)) {
            if (fisherMatrix.getData() == 0) 
                fisherMatrix = ndarray::allocate(modelMatrix.cols(), modelMatrix.cols());
            ndarray::EigenView<double,2,2> fisherMatrixE(fisherMatrix);
            fisherMatrixE.part<Eigen::SelfAdjoint>() = modelMatrix.transpose() * modelMatrix;
            fisherMatrixE.part<Eigen::SelfAdjoint>() += _priorMatrix.part<Eigen::SelfAdjoint>();
            Bit<COEFFICIENT_FISHER_MATRIX>::set(status);
        }
        _llt.compute(ndarray::viewAsEigen(fisherMatrix));
        fisherFactor = _llt.matrixL();
    }

    virtual void computeFisherMatrix(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double,2,2> fisherMatrix,
        ndarray::Array<double,2,2> & fisherFactor,
        int & status
    ) {
        if (Bit<COEFFICIENT_FISHER_FACTOR>::test(status)) {
            fisherMatrix.part<Eigen::SelfAdjoint>() 
                = ndarray::viewAsEigen(fisherFactor).part<Eigen::LowerTriangular>()
                * ndarray::viewAsTransposedEigen(fisherFactor).part<Eigen::UpperTriangular>();
        } else {
            fisherMatrix.part<Eigen::SelfAdjoint>() = modelMatrix.transpose() * modelMatrix;
            fisherMatrix.part<Eigen::SelfAdjoint>() += _priorMatrix.part<Eigen::SelfAdjoint>();
        }
    }

    virtual void setPrior(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma) {
        _llt.compute(sigma);
        _priorMatrix = MatrixRM::Identity(sigma.rows(), sigma.cols());
        _llt.solveInPlace(_priorMatrix);
        _priorVector = _priorMatrix * mu;
    }

    virtual void setPrior(Eigen::VectorXd const & mu) {
        _priorVector = _priorMatrix * mu;
    }

    explicit CholeskySolver(int size) : 
        _priorMatrix(MatrixRM::Zero(size, size)),
        _priorVector(Eigen::VectorXd::Zero(size)),
        _rhs(Eigen::VectorXd::Zero(size))
    {}

private:
    LLT _llt;
    MatrixRM _priorMatrix;
    Eigen::VectorXd _priorVector;
    Eigen::VectorXd _rhs;
};

class Evaluation::EigenSolver : public Evaluation::LinearSolver {
public:
    typedef Eigen::SelfAdjointEigenSolver<MatrixRM> Factorization;

    virtual void solve(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double const,1,1> const & data,
        ndarray::EigenView<double,1,1> coefficients,
        ndarray::Array<double,2,2> & fisherMatrix,
        ndarray::Array<double,2,2> & fisherFactor,
        int & status
    ) {
        if (!Bit<COEFFICIENT_FISHER_MATRIX>::test(status)) {
            if (fisherMatrix.getData() == 0)
                fisherMatrix = ndarray::allocate(modelMatrix.cols(), modelMatrix.cols());
            computeFisherMatrix(modelMatrix, ndarray::viewAsEigen(fisherMatrix), fisherFactor, status);
            Bit<COEFFICIENT_FISHER_MATRIX>::set(status);
        }
        _rhs = (modelMatrix.transpose() * data).lazy();
        _rhs += _priorVector;
        _factorization.compute(ndarray::viewAsEigen(fisherMatrix), true);
        Eigen::VectorXd values = _factorization.eigenvalues();
        double rcond = values.maxCoeff() * std::sqrt(std::numeric_limits<double>::epsilon());
        values = (values.cwise() < rcond).select(
            Eigen::VectorXd::Zero(values.size()), values.cwise().inverse()
        );
        coefficients = _factorization.eigenvalues() * values.asDiagonal()
            * _factorization.eigenvalues().transpose() * _rhs;        
    }

    virtual void computeFisherFactor(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double,2,2> fisherFactor,
        ndarray::Array<double,2,2> & fisherMatrix,
        int & status
    ) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Robust solvers do not compute the Fisher matrix Cholesky factorization."
        );
    }

    virtual void computeFisherMatrix(
        ndarray::EigenView<double const,2,2> const & modelMatrix,
        ndarray::EigenView<double,2,2> fisherMatrix,
        ndarray::Array<double,2,2> & fisherFactor,
        int & status
    ) {
        fisherMatrix.part<Eigen::SelfAdjoint>() = modelMatrix.transpose() * modelMatrix;
        fisherMatrix.part<Eigen::SelfAdjoint>() += _priorMatrix.part<Eigen::SelfAdjoint>();
    }

    virtual void setPrior(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma) {
        _factorization.compute(sigma);
        _priorMatrix = _factorization.eigenvectors() 
            * _factorization.eigenvalues().cwise().inverse().asDiagonal()
            * _factorization.eigenvectors().transpose();
        _priorVector = _priorMatrix * mu;
    }

    virtual void setPrior(Eigen::VectorXd const & mu) {
        _priorVector = _priorMatrix * mu;
    }

    explicit EigenSolver(int size) : 
        _priorMatrix(MatrixRM::Zero(size, size)),
        _priorVector(Eigen::VectorXd::Zero(size)),
        _rhs(Eigen::VectorXd::Zero(size))
    {}

private:

    Factorization _factorization;
    MatrixRM _priorMatrix;
    Eigen::VectorXd _priorVector;
    Eigen::VectorXd _rhs;
};

Evaluation::Evaluation(BaseEvaluator::Ptr const & evaluator, bool robustSolver) : 
    _status(0), _evaluator(evaluator), _parameters(ndarray::allocate(_evaluator->getParameterSize()))
{
    _evaluator->writeInitialParameters(_parameters);
    _solver.reset(new CholeskySolver(evaluator->getCoefficientSize()));
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator, BaseDistribution const & prior, bool robustSolver
) : 
    _status(0), _evaluator(evaluator), _prior(prior.clone()),
    _parameters(ndarray::allocate(_evaluator->getParameterSize()))
{
    _evaluator->writeInitialParameters(_parameters);
    initialize();
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator,
    lsst::ndarray::Array<double const,1,1> const & parameters,
    bool robustSolver
) : 
    _status(0), _evaluator(evaluator), _parameters(ndarray::copy(parameters))
{
    initialize();
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator,
    lsst::ndarray::Array<double const,1,1> const & parameters,
    BaseDistribution const & prior,
    bool robustSolver
) : 
    _status(0), _evaluator(evaluator), _prior(prior.clone()),
    _parameters(ndarray::copy(parameters))
{
    initialize();
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator,
    Eigen::VectorXd const & parameters,
    bool robustSolver
) : 
    _status(0), _evaluator(evaluator), _parameters(ndarray::allocate(parameters.size()))
{
    ndarray::viewAsEigen(_parameters) = parameters;
    initialize();
}

Evaluation::Evaluation(
    BaseEvaluator::Ptr const & evaluator,
    Eigen::VectorXd const & parameters,
    BaseDistribution const & prior,
    bool robustSolver
) : 
    _status(0), _evaluator(evaluator), _prior(prior.clone()),
    _parameters(ndarray::allocate(parameters.size()))
{
    ndarray::viewAsEigen(_parameters) = parameters;
    initialize();
}

void Evaluation::update(lsst::ndarray::Array<double const,1,1> const & parameters) {
    assert(parameters.getSize<0>() == _parameters.getSize<0>());
    _parameters.deep() = parameters;
    _status = 0;
    updateNestedPrior();
}
void Evaluation::update(Eigen::VectorXd const & parameters) {
    assert(parameters.size() == _parameters.getSize<0>());
    ndarray::viewAsEigen(_parameters) = parameters;
    _status = 0;
    updateNestedPrior();
}
    
void Evaluation::update(
    lsst::ndarray::Array<double const,1,1> const & parameters, 
    lsst::ndarray::Array<double const,1,1> const & coefficients
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

void Evaluation::setCoefficients(lsst::ndarray::Array<double const,1,1> const & coefficients) {
    assert(coefficients.getSize<0>() == _coefficients.getSize<0>());
    _coefficients.deep() = coefficients;
    _status &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::set(_status);
}

void Evaluation::setCoefficients(Eigen::VectorXd const & coefficients) {
    assert(coefficients.size() == _coefficients.getSize<0>());
    ndarray::viewAsEigen(_coefficients) = coefficients;
    _status &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::set(_status);
}

void Evaluation::solveCoefficients() {
    _status &= ~coefficient_dependencies;
    Bit<COEFFICIENTS>::reset(_status);
    ensureCoefficients();
}

void Evaluation::ensureModelMatrix() const {
    if (Bit<MODEL_MATRIX>::test(_status)) return;
    if (_modelMatrix.getData() == 0) {
        _modelMatrix = ndarray::allocate(_evaluator->getDataSize(), _evaluator->getCoefficientSize());
    }
    _evaluator->_evaluateModelMatrix(_modelMatrix, _parameters);
    Bit<MODEL_MATRIX>::set(_status);
}

void Evaluation::ensureModelMatrixDerivative() const {
    if (Bit<MODEL_MATRIX_DERIVATIVE>::test(_status)) return;
    ensureModelMatrix();
    if (_modelMatrixDerivative.getData() == 0) {
        _modelMatrixDerivative = ndarray::allocate(
            _evaluator->getParameterSize(), _evaluator->getDataSize(), _evaluator->getCoefficientSize()
        );
    }
    _evaluator->_evaluateModelMatrixDerivative(_modelMatrixDerivative, _modelMatrix, _parameters);
    Bit<MODEL_MATRIX_DERIVATIVE>::set(_status);
}

void Evaluation::ensureCoefficients() const {
    if (Bit<COEFFICIENTS>::test(_status)) return;
    ensureModelMatrix();
    if (_coefficients.getData() == 0) {
        _coefficients = ndarray::allocate(_evaluator->getCoefficientSize());
    }
    _solver->solve(
        ndarray::EigenView<double const,2,2>(_modelMatrix),
        ndarray::EigenView<double const,1,1>(_evaluator->getDataVector()),
        ndarray::EigenView<double,1,1>(_coefficients),
        _coefficientFisherMatrix,
        _coefficientFisherFactor,
        _status
    );
    Bit<COEFFICIENTS>::set(_status);
}

void Evaluation::ensureResiduals() const {
    if (Bit<RESIDUALS>::test(_status)) return;
    ensureCoefficients();
    ensureModelMatrix();
    if (_residuals.getData() == 0) {
        _residuals = ndarray::allocate(_evaluator->getDataSize());
    }
    ndarray::viewAsEigen(_residuals) 
        = ndarray::viewAsEigen(_modelMatrix) * ndarray::viewAsEigen(_coefficients);
    ndarray::viewAsEigen(_residuals) -= ndarray::viewAsEigen(_evaluator->getDataVector());
    Bit<RESIDUALS>::set(_status);
}

void Evaluation::ensureResidualsJacobian() const {
    if (Bit<RESIDUALS_JACOBIAN>::test(_status)) return;
    ensureCoefficients();
    ensureModelMatrixDerivative();
    if (_residualsJacobian.getData() == 0) {
        _residualsJacobian = ndarray::allocate(
            _evaluator->getDataSize(), _evaluator->getParameterSize()
        );
    }
    ndarray::EigenView<double,1,1> coeffVec(_coefficients);
    ndarray::EigenView<double,2,2> jac(_residualsJacobian);
    for (int n = 0; n < jac.cols(); ++n) {
        jac.col(n) = ndarray::viewAsEigen(_modelMatrixDerivative[n]) * coeffVec;
    }
    Bit<RESIDUALS_JACOBIAN>::set(_status);
}

void Evaluation::ensureObjectiveValue() const {
    if (Bit<OBJECTIVE_VALUE>::test(_status)) return;
    ensureResiduals();
    _objectiveValue = 0.5 * ndarray::viewAsEigen(_residuals).squaredNorm();
    if (_prior) {
        _objectiveValue -= std::log(_prior->evaluate(_parameters));
        if (_nestedPrior) {
            _objectiveValue -= std::log(_nestedPrior->getNormalization());
        }
    }
    Bit<OBJECTIVE_VALUE>::set(_status);
}

void Evaluation::ensureCoefficientFisherMatrix() const {
    if (Bit<COEFFICIENT_FISHER_MATRIX>::test(_status)) return;
    ensureModelMatrix();
    if (_coefficientFisherMatrix.getData() == 0) {
        _coefficientFisherMatrix = ndarray::allocate(
            _evaluator->getCoefficientSize(), _evaluator->getCoefficientSize()
        );
    }
    _solver->computeFisherMatrix(
        ndarray::EigenView<double const,2,2>(_modelMatrix),
        ndarray::EigenView<double,2,2>(_coefficientFisherMatrix),
        _coefficientFisherFactor,
        _status
    );
    Bit<COEFFICIENT_FISHER_MATRIX>::set(_status);
}

void Evaluation::ensureCoefficientFisherFactor() const {
    if (Bit<COEFFICIENT_FISHER_FACTOR>::test(_status)) return;
    ensureModelMatrix();
    if (_coefficientFisherFactor.getData() == 0) {
        _coefficientFisherFactor = ndarray::allocate(
            _evaluator->getCoefficientSize(), _evaluator->getCoefficientSize()
        );
    }
    _solver->computeFisherFactor(
        ndarray::EigenView<double const,2,2>(_modelMatrix),
        ndarray::EigenView<double,2,2>(_coefficientFisherFactor),
        _coefficientFisherMatrix,
        _status
    );
    Bit<COEFFICIENT_FISHER_FACTOR>::set(_status);
}

void Evaluation::initialize() {
    if (_evaluator->getParameterSize() != _parameters.getSize<0>()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            (boost::format("Evaluator parameter size (%d) does not match parameter vector size (%d).")
             % _evaluator->getParameterSize() % _parameters.getSize<0>()).str()
        );
    }
    if (_prior) {
        if (_evaluator->getParameterSize() != _prior->getDimensionality()) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LengthErrorException,
                (boost::format("Evaluator parameter size (%d) does not match prior dimensionality (%d).")
                 % _evaluator->getParameterSize() % _prior->getDimensionality()).str()
            );
        }
        BaseDistribution::Ptr nestedAsBase = _prior->evaluateNested(_parameters);
        if (nestedAsBase) {
            _nestedPrior = boost::dynamic_pointer_cast<GaussianDistribution>(nestedAsBase);
            if (!_nestedPrior) {
                throw LSST_EXCEPT(
                    lsst::pex::exceptions::InvalidParameterException,
                    "Evaluation only supports prior distributions with no nesting or "
                "a nested GaussianDistribution."
                );
            }
            if (_evaluator->getCoefficientSize() != _nestedPrior->getDimensionality()) {
                throw LSST_EXCEPT(
                    lsst::pex::exceptions::LengthErrorException,
                    (boost::format("Evaluator coefficient size (%d) does not match prior "
                                   "nested dimensionality (%d).")
                     % _evaluator->getCoefficientSize() % _nestedPrior->getDimensionality()).str()
                );
            }
            _solver.reset(new CholeskySolver(_evaluator->getCoefficientSize()));
            _solver->setPrior(_nestedPrior->getMu(), _nestedPrior->getSigma());
        }
    }
    if (!_solver) {
        // TODO: someday this could be a QR solver here.
        _solver.reset(new CholeskySolver(_evaluator->getCoefficientSize()));
    }
}

void Evaluation::updateNestedPrior() {
    if (_nestedPrior) {
        int depFlags = _prior->getNestedDependency();
        if (depFlags) {
            _prior->updateNested(*_nestedPrior, _parameters);
            if (depFlags & GaussianDistribution::SIGMA_DEPENDENT) {
                _solver->setPrior(_nestedPrior->getMu(), _nestedPrior->getSigma());
            } else {
                _solver->setPrior(_nestedPrior->getMu());
            }
        }
    }
}

Evaluation::~Evaluation() {}


}}} // namespace lsst::meas::multifit
