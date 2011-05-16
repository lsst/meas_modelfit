#include "lsst/meas/multifit/qp.h"
#include "lsst/pex/exceptions.h"
#include "boost/make_shared.hpp"
#include "QuadProg++.hh"
#undef solve
namespace lsst { namespace meas { namespace multifit {

namespace {

QuadProgPP::Matrix<double> makeMatrix(lsst::ndarray::Array<double const,2,0> const & array) {
    QuadProgPP::Matrix<double> r(array.getSize<0>(), array.getSize<1>());
    lsst::ndarray::Array<double const,2,0>::Iterator i = array.begin();
    for (int n = 0; i != array.end(); ++n, ++i) {
        std::copy(i->begin(), i->end(), r[n]);
    }
    return r;
}

QuadProgPP::Vector<double> makeVector(lsst::ndarray::Array<double const,1,0> const & array) {
    QuadProgPP::Vector<double> r(array.getSize<0>());
    lsst::ndarray::Array<double const,1,0>::Iterator i = array.begin();
    for (int n = 0; i != array.end(); ++n, ++i) r[n] = *i;
    return r;
}

} // anonymous

struct QPSolver::Data {
    QuadProgPP::Matrix<double> g;
    QuadProgPP::Vector<double> c;
    QuadProgPP::Matrix<double> ae;
    QuadProgPP::Vector<double> be;
    QuadProgPP::Matrix<double> ai;
    QuadProgPP::Vector<double> bi;
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
    _data->g = makeMatrix(g);
    _data->c = makeVector(c);
    _data->ae.resize(n, 0);
    _data->ai.resize(n, 0);
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
    _data->ae = makeMatrix(a.transpose());
    _data->be = makeVector(b); 
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
    _data->ai = makeMatrix(a.transpose());
    _data->bi = makeVector(b); 
    _data->bi *= -1.0;
    return *this;
}

double QPSolver::solve(lsst::ndarray::Array<double,1,0> const & x) const {
    QuadProgPP::Vector<double> qp_x;
    double cost = QuadProgPP::solve_quadprog(
        _data->g, _data->c, _data->ae, _data->be, _data->ai, _data->bi, qp_x
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
