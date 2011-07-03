#include "lsst/meas/multifit/qp.h"
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
    if (g.getSize<0>() != n) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of QP matrix rows does not match vector size."
        );
    }
    if (g.getSize<1>() != n) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of QP matrix columns does not match vector size."
        );
    }
    _data->g = ndarray::viewAsEigen(g);
    _data->c = ndarray::viewAsEigen(c);
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
    if (a.getSize<0>() != b.getSize<0>()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of rows of constraint matrix does not match number of rows of vector."
        );
    }
    _data->ae = ndarray::viewAsTransposedEigen(a);
    _data->be = ndarray::viewAsEigen(b); 
    _data->be *= -1.0;
    return *this;
}

QPSolver & QPSolver::inequality(
    lsst::ndarray::Array<double const,2,0> const & a,
    lsst::ndarray::Array<double const,1,0> const & b
) {
    if (a.getSize<0>() != b.getSize<0>()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of rows of constraint matrix does not match number of rows of vector."
        );
    }
    _data->ai = ndarray::viewAsTransposedEigen(a);
    _data->bi = ndarray::viewAsEigen(b); 
    _data->bi *= -1.0;
    return *this;
}

double QPSolver::solve(lsst::ndarray::Array<double,1,0> const & x) const {
    Eigen::VectorXd qp_x;
    double cost = QuadProgPP::solve_quadprog(
        _data->g, _data->c, _data->ae, _data->be, _data->ai, _data->bi, qp_x,
	_data->maxIterations
    );
    if (x.getSize<0>() != int(qp_x.size())) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Incorrect size for QP solution vector."
        );
    }
    for (int n = 0; n < x.getSize<0>(); ++n) {
        x[n] = qp_x[n];
    }
    return cost;
}

}}} // namespace lsst::meas::multifit
