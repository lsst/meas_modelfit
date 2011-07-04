#include "lsst/meas/multifit/qp.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/pex/exceptions.h"
#include "boost/make_shared.hpp"
#include "QuadProg++.hh"
#include "lsst/ndarray/eigen.h"

#undef solve
namespace lsst { namespace meas { namespace multifit {

struct QPSolver::Data {
    int maxIterations;
    Eigen::MatrixXd g;
    Eigen::VectorXd c;
    Eigen::MatrixXd ae;
    Eigen::VectorXd be;
    Eigen::MatrixXd ai;
    Eigen::VectorXd bi;
};

QPSolver::QPSolver(
    lsst::ndarray::Array<double const,2,0> const & g,
    lsst::ndarray::Array<double const,1,0> const & c
) : _data(boost::make_shared<Data>())
{
    int n = c.getSize<0>();
    detail::checkSize(
        g.getSize<0>(), n, 
        "Number of QP matrix rows (%d) does not match vector size (%d)."
    );
    detail::checkSize(
        g.getSize<1>(), n, 
        "Number of QP matrix columns (%d) does not match vector size (%d)."
    );
    _data->g = ndarray::viewAsEigen(g);
    _data->c = ndarray::viewAsEigen(c);
    _data->maxIterations = 100;
}

QPSolver::QPSolver(
    Eigen::MatrixXd const & g,
    Eigen::VectorXd const & c
) : _data(boost::make_shared<Data>())
{
    int n = c.size();
    detail::checkSize(
        g.rows(), n, 
        "Number of QP matrix rows (%d) does not match vector size (%d)."
    );
    detail::checkSize(
        g.cols(), n, 
        "Number of QP matrix columns (%d) does not match vector size (%d)."
    );
    _data->g = g;
    _data->c = c;
    _data->maxIterations = 100;
}

QPSolver & QPSolver::maxIterations(int n) { 
    _data->maxIterations = n;
    return *this;
}

QPSolver & QPSolver::equality(
    lsst::ndarray::Array<double const,2,0> const & a,
    lsst::ndarray::Array<double const,1,0> const & b
) { 
    detail::checkSize(
        a.getSize<0>(), b.getSize<0>(),
        "Number of rows of constraint matrix (%d) does not match size of vector (%d)."
    );
    detail::checkSize(
        a.getSize<1>(), _data->c.size(),
        "Number of columns of constraint matrix (%d) does not match size of problem (%d)."
    );
    _data->ae = ndarray::viewAsTransposedEigen(a);
    _data->be = ndarray::viewAsEigen(b); 
    _data->be *= -1.0;
    return *this;
}

QPSolver & QPSolver::equality(
    Eigen::MatrixXd const & a,
    Eigen::VectorXd const & b
) { 
    detail::checkSize(
        a.rows(), b.size(),
        "Number of rows of constraint matrix (%d) does not match size of vector (%d)."
    );
    detail::checkSize(
        a.cols(), _data->c.size(),
        "Number of columns of constraint matrix (%d) does not match size of problem (%d)."
    );
    _data->ae = a.transpose();
    _data->be = b;
    _data->be *= -1.0;
    return *this;
}

QPSolver & QPSolver::inequality(
    lsst::ndarray::Array<double const,2,0> const & a,
    lsst::ndarray::Array<double const,1,0> const & b
) {
    detail::checkSize(
        a.getSize<0>(), b.getSize<0>(),
        "Number of rows of constraint matrix (%d) does not match size of vector (%d)."
    );
    detail::checkSize(
        a.getSize<1>(), _data->c.size(),
        "Number of columns of constraint matrix (%d) does not match size of problem (%d)."
    );
    _data->ai = ndarray::viewAsTransposedEigen(a);
    _data->bi = ndarray::viewAsEigen(b); 
    _data->bi *= -1.0;
    return *this;
}

QPSolver & QPSolver::inequality(
    Eigen::MatrixXd const & a,
    Eigen::VectorXd const & b
) {
    detail::checkSize(
        a.rows(), b.size(),
        "Number of rows of constraint matrix (%d) does not match size of vector (%d)."
    );
    detail::checkSize(
        a.cols(), _data->c.size(),
        "Number of columns of constraint matrix (%d) does not match size of problem (%d)."
    );
    _data->ai = a.transpose();
    _data->bi = b;
    _data->bi *= -1.0;
    return *this;
}

double QPSolver::solve(lsst::ndarray::Array<double,1,0> const & x) const {
    Eigen::VectorXd qp_x;
    double cost = QuadProgPP::solve_quadprog(
        _data->g, _data->c, _data->ae, _data->be, _data->ai, _data->bi, qp_x,
	_data->maxIterations
    );
    detail::checkSize(
        x.getSize<0>(), qp_x.size(),
        "Size of solution vector (%d) is incorrect (%d)."
    );
    ndarray::viewAsEigen(x) = qp_x;
    return cost;
}

double QPSolver::solve(Eigen::VectorXd & x) const {
    double cost = QuadProgPP::solve_quadprog(
        _data->g, _data->c, _data->ae, _data->be, _data->ai, _data->bi, x,
	_data->maxIterations
    );
    return cost;
}

}}} // namespace lsst::meas::multifit
