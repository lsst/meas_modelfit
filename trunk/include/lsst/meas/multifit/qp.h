// -*- LSST-C++ -*-
#ifndef LSST_MEAS_MULTIFIT_qp_h_INCLUDED
#define LSST_MEAS_MULTIFIT_qp_h_INCLUDED

#include "lsst/ndarray.h"
#include "boost/shared_ptr.hpp"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Solve a quadratic programming problem using QuadProg++.
 *
 *  @f$ \min_x q(x) = \frac{1}{2} x^T G x + x^T c @f$
 *  subject to @f$ A_e x = b_e @f$ and @f$ A_i x = b_i @f$.
 */
class QPSolver {
public:

    /**
     *  @brief Initialize the solver.
     *
     *  @param[in]      g   The @f$N \times N@f$ matrix @f$G@f$.
     *  @param[in]      c   The @f$N@f$-element vector @f$c@f$.
     */
    QPSolver(
        lsst::ndarray::Array<double const,2,0> const & g,
        lsst::ndarray::Array<double const,1,0> const & c
    );

    /**
     *  @brief Set the equality constraint for the solver.
     *
     *  @param[in]      a   The @f$K_e \times N@f$ matrix @f$A_e@f$.
     *  @param[in]      b   The @f$K_e@f$-element vector @f$b_e@f$.
     *
     *  Multiple calls will reset the equality constraint, not append to it.
     */
    QPSolver & equality(
        lsst::ndarray::Array<double const,2,0> const & a,
        lsst::ndarray::Array<double const,1,0> const & b
    );

    /**
     *  @brief Set the inequality constraint for the solver.
     *
     *  @param[in]      a   The @f$K_i \times N@f$ matrix @f$A_i@f$.
     *  @param[in]      b   The @f$K_i@f$-element vector @f$b_i@f$.
     *
     *  Multiple calls will reset the inequality constraint, not append to it.
     */
    QPSolver & inequality(
        lsst::ndarray::Array<double const,2,0> const & a,
        lsst::ndarray::Array<double const,1,0> const & b
    );

    /**
     *  @brief Solve the quadratic programming problem.
     *
     *  @param[in,out]  x   The solution vector @f$x@f$.  Must be allocated to the correct size,
     *                      but the values of the array elements will be ignored.
     *
     *  @return the value of the cost function @f$q(x)@f$ at the solution, or infinity if the
     *          problem is infeasible.
     */
    double solve(lsst::ndarray::Array<double,1,0> const & x) const;

private:
    struct Data;
    boost::shared_ptr<Data> _data;
};

}}} // namespace lsst::meas::multifit

#endif // LSST_MEAS_MULTIFIT_qp_h_INCLUDED
