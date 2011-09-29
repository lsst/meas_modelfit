// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 *
 * Copyright 2011 LSST Corporation.
 * Copyright 2007-2009 Luca Di Gaspero.
 * Copyright 2009 Eric Moyer.  
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 *
 * Much of the code in this file originated as part of QuadProg++: 
 * a C++ library implementing the algorithm of Goldfarb and Idnani for
 * the solution of a (convex) Quadratic Programming problem by means
 * of an active-set dual method, by Luca Di Gaspero (LGPL v3).  It was
 * heavily modified and adapted for use in the LSST software framework
 * by Jim Bosch in 2011.
 */

#include "lsst/meas/multifit/qp.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/pex/exceptions.h"
#include "boost/make_shared.hpp"
#include "lsst/ndarray/eigen.h"
#include "Eigen/Cholesky"

#undef solve
namespace lsst { namespace meas { namespace multifit {

struct QPSolver::Data {
    int maxIterations;
    bool isDiagonal;
    Eigen::MatrixXd g;
    Eigen::VectorXd c;
    Eigen::MatrixXd ae;
    Eigen::VectorXd be;
    Eigen::MatrixXd ai;
    Eigen::VectorXd bi;
};

class QPSolver::Impl {
public:
    Data const & input;
    Eigen::VectorXd & x;

    int const n; // size of problem
    int const p; // number of equality constraints
    int const m; // numer of inequality constraints

    static double const INF;
    static double const EPS;

    Impl(Data const & input_, Eigen::VectorXd & x_);

    double solve();

    Eigen::VectorXd d;
    Eigen::VectorXd z;
    Eigen::VectorXd r;
    Eigen::VectorXd u;
    Eigen::VectorXi A;
    Eigen::MatrixXd J;
    Eigen::MatrixXd R;

    // In-place lower Cholesky factorization.
    static void factorCholesky(Eigen::MatrixXd & matrix, bool isDiagonal);

    // Solve @f$ L y = b @f$
    static void subForward(
        Eigen::MatrixXd const & L, Eigen::VectorXd & y, Eigen::VectorXd const & b, bool isDiagonal
    );

    // Solve @f$ L^T x = y @f$
    static void subBackward(
        Eigen::MatrixXd const & L, Eigen::VectorXd & x, Eigen::VectorXd const & y, bool isDiagonal
    );

    // Solve @f$ L L^T x = b @f$
    static void solveCholesky(
        Eigen::MatrixXd const & L, Eigen::VectorXd & x, Eigen::VectorXd const & b, bool isDiagonal
    );

    // Set @f$ J = L^{-T} @f$
    static void initializeJ(Eigen::MatrixXd const & L, Eigen::MatrixXd & J, bool isDiagonal);

    static inline double distance(double a, double b) {
        register double a1, b1, t;
        a1 = std::abs(a);
        b1 = std::abs(b);
        if (a1 > b1) {
            t = (b1 / a1);
            return a1 * ::std::sqrt(1.0 + t * t);
        } else if (b1 > a1) {
            t = (a1 / b1);
            return b1 * ::std::sqrt(1.0 + t * t);
        }
        return a1 * ::std::sqrt(2.0);
    }

    void update(int iq);
    bool add_constraint(int& iq, double& rnorm);
    void delete_constraint(int& iq, int l);
};

double const QPSolver::Impl::INF = std::numeric_limits<double>::infinity();
double const QPSolver::Impl::EPS = std::numeric_limits<double>::epsilon();

QPSolver::Impl::Impl(Data const & input_, Eigen::VectorXd & x_) :
    input(input_), x(x_), n(input.c.size()), p(input.be.size()), m(input.bi.size()),
    d(Eigen::VectorXd::Zero(n)),
    z(n),
    r(m + p),
    u(m + p),
    A(m + p),
    J(Eigen::MatrixXd::Identity(n, n)),
    R(Eigen::MatrixXd::Zero(n, n))
{
    x.resize(n);
}

double QPSolver::Impl::solve() {
    Eigen::MatrixXd G(input.g);
    factorCholesky(G, input.isDiagonal);
    double R_norm = 1.0;
    initializeJ(G, J, input.isDiagonal);
    double c = input.g.trace() * J.trace(); // c is an estimate of cond(G)

    // Solve for unconstrained minimizer of 0.5 * x^T G x + c^T x
    solveCholesky(G, x, input.c, input.isDiagonal);
    x *= -1;
    double objective = 0.5 * input.c.dot(x);

    // Add equality constraints to working set.
    int iq = 0;
    Eigen::VectorXd x_old(n), u_old(m + p);
    Eigen::VectorXi A_old(m + p), iai(m + p);
    for (int i = 0; i < p; ++i) {
        d = J.transpose() * input.ae.col(i);
        update(iq);
        // Compute full step length t2: i.e., the minimum step in primal space s.t. the contraint 
        // becomes feasible
        double t2 = 0.0;
        if (std::abs(z.squaredNorm()) > EPS) // i.e. z != 0
            t2 = (-input.ae.col(i).dot(x) - input.be[i]) / z.dot(input.ae.col(i));
        x += t2 * z;
    
        /* set u = u+ */
        u[iq] = t2;
        if (iq > 0)
            u.segment(0, iq) -= t2 * r.segment(0, iq);
    
        /* compute the new solution value */
        objective += 0.5 * (t2 * t2) * z.dot(input.ae.col(i));
        A[i] = -i - 1;
    
        if (!add_constraint(iq, R_norm)) {	  
            // Equality constraints are linearly dependent
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Constraints are linearly dependent"
            );
            return objective;
        }
    }

    for (int i =0; i < m; ++i) iai[i] = i;

    int iterations = 0;
    Eigen::VectorXi iaexcl(m + p);
    Eigen::VectorXd s(m + p);
    int ip = 0;
    while (true) {
        // step 1: choose a violated constraint
        for (int i = p; i < iq; i++) {
            ip = A[i];
            iai[ip] = -1;
        }
	
        // compute s[x] = ai^T * x + bi for all elements of K \ A
        double ss = 0.0;
        double psi = 0.0; // this value will contain the sum of all infeasibilities
        ip = 0;    // ip will be the index of the chosen violated constraint 
        s = input.ai.transpose() * x + input.bi;
        for (int i = 0; i < m; ++i) {
            iaexcl[i] = true;
            psi += std::min(0.0, s[i]);
        }
  
        if (std::abs(psi) <= m * EPS * c * 100.0) {
            // numerically there are not infeasibilities anymore
            return objective;
        }
 
        // save old values for u and A
        if (iq > 0) {
            u_old.segment(0, iq) = u.segment(0, iq);
            A_old.segment(0, iq) = A.segment(0, iq);
        }
        // and for x
        x_old = x;

        while (true) { // Step 2: check for feasibility and determine a new S-pair
            for (int i = 0; i < m; ++i) {
                if (s[i] < ss && iai[i] != -1 && iaexcl[i]) {
                    ss = s[i];
                    ip = i;
                }
            }
            if (ss >= 0.0) {
                return objective;
            }
    
            // set u = [u 0]^T
            u[iq] = 0.0;
            // add ip to the active set A
            A[iq] = ip;

            double t, t1, t2;
            bool repeat1 = false;
            while (true) { // Step 2a: determine step direction
                if (++iterations > input.maxIterations) {
                    throw LSST_EXCEPT(
                        lsst::pex::exceptions::RuntimeErrorException,
                        "Exceeded maximum number of QP iterations."
                    );
                }
                // compute z = H np: the step direction in the primal space (through J, see the paper)
                // and  N* np (if q > 0): the negative of the step direction in the dual space
                d = J.transpose() * input.ai.col(ip);
                update(iq);
                
                // Step 2b: compute step length
                int l = 0;
                // Compute t1: partial step length (maximum step in dual space without 
                // violating dual feasibility
                t1 = INF;
                /* find the index l s.t. it reaches the minimum of u+[x] / r */
                for (int k = p; k < iq; k++) {
                    if (r[k] > 0.0) {
                        if (u[k] / r[k] < t1) {
                            t1 = u[k] / r[k];
                            l = A[k];
                        }
                    }
                }
                // Compute t2: full step length (minimum step in primal space such that the
                // constraint ip becomes feasible */
                t2 = INF;
                if (std::fabs(z.squaredNorm()) > EPS) // i.e. z != 0
                    t2 = -s[ip] / z.dot(input.ai.col(ip));
  
                // the step is chosen as the minimum of t1 and t2
                t = std::min(t1, t2);

                // Step 2c: determine new S-pair and take step:
  
                // case (i): no step in primal or dual space
                if (t >= INF) {
                    /* QPP is infeasible */
                    return INF;
                }

                // case (ii): step in dual space
                if (t2 >= INF) {
                    /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
                    if (iq > 0)
                        u.segment(0, iq) -= t * r.segment(0, iq);
                    u[iq] += t;
                    iai[l] = l;
                    delete_constraint(iq, l);
                    // repeat Step 2a
                    continue;
                }

                /* case (iii): step in primal and dual space */
  
                /* set x = x + t * z */
                x = x + t * z;
                /* update the solution value */
                objective += t * z.dot(input.ai.col(ip)) * (0.5 * t + u[iq]);
                /* u = u + t * [-r 1] */
                if (iq > 0)
                    u.segment(0, iq) -= t * r.segment(0, iq);
                u[iq] += t;

                if (std::fabs(t - t2) < EPS) {
                    // full step has taken
                    // add constraint ip to the active set
                    if (!add_constraint(iq, R_norm)) {
                        iaexcl[ip] = false;
                        delete_constraint(iq, ip);
                        for (int i = 0; i < m; i++)
                            iai[i] = i;
                        for (int i = p; i < iq; i++) {
                            A[i] = A_old[i];
                            u[i] = u_old[i];
                            iai[A[i]] = -1;
                        }
                        x = x_old;
                        // repeat Step 2
                        break;
                    } else {
                        iai[ip] = -1;
                    }
                    repeat1 = true;
                    break;
                }
                
                // drop constraint l
                iai[l] = l;
                delete_constraint(iq, l);

                // update s[ip] = CI * x + ci0
                s[ip] = input.ai.col(ip).dot(x) + input.bi[ip];

            } // while for Step 2a
            if (repeat1) break;

        } // while for Step 2

    } // while for Step 1
    
    return INF;
}

void QPSolver::Impl::factorCholesky(Eigen::MatrixXd & matrix, bool isDiagonal) {
    if (isDiagonal) {
        matrix.diagonal().array() = matrix.diagonal().array().sqrt();
    } else {
        Eigen::LLT<Eigen::MatrixXd> cholesky(matrix);
        if (!cholesky.isPositiveDefinite()) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::RuntimeErrorException, 
                "Singular matrix in QP Cholesky decomposition."
            );
        }
        matrix = cholesky.matrixL();
    }
}

void QPSolver::Impl::subForward(
    Eigen::MatrixXd const & L, Eigen::VectorXd & y, Eigen::VectorXd const & b, bool isDiagonal
) {
    if (isDiagonal) {
        y.array() = b.array() / L.diagonal().array();
    } else {
        y = L.triangularView<Eigen::Lower>().solve(b);
    }
}

void QPSolver::Impl::subBackward(
    Eigen::MatrixXd const & L, Eigen::VectorXd & x, Eigen::VectorXd const & y, bool isDiagonal
) {
    if (isDiagonal) {
        x.array() = y.array() / L.diagonal().array();
    } else {
        x = L.triangularView<Eigen::Lower>().transpose().solveTriangular(y);
    }
}

void QPSolver::Impl::solveCholesky(
    Eigen::MatrixXd const & L, Eigen::VectorXd & x, Eigen::VectorXd const & b, bool isDiagonal
) {
    if (isDiagonal) {
        x.array() = b.array() / L.diagonal().array().square();
    } else {
        x = L.triangularView<Eigen::Lower>().solve(b);
        L.triangularView<Eigen::Lower>().transpose().solveInPlace(x);
    }
}

void QPSolver::Impl::initializeJ(Eigen::MatrixXd const & L, Eigen::MatrixXd & J, bool isDiagonal) {
    J.setIdentity();
    if (isDiagonal) {
        J.diagonal().array() = L.diagonal().array().inverse();
    } else {
        L.triangularView<Eigen::Lower>().transpose().solveInPlace(J);
    }
}

void QPSolver::Impl::update(int iq) {
    /* set z = H * d */
    if (iq < n)
        z = J.block(0, iq, n, n - iq) * d.segment(iq, n - iq); 
    else
	z.setZero();

    /* set r = R^-1 d */
    for (int i = iq - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < iq; ++j)
            sum += R(i,j) * r[j];
        r[i] = (d[i] - sum) / R(i,i);
    }
}

bool QPSolver::Impl::add_constraint(int& iq, double& R_norm) { 
    double cc, ss, h, t1, t2, xny;
	
    /* we have to find the Givens rotation which will reduce the element
       d[j] to zero.
       if it is already zero we don't have to do anything, except of
       decreasing j */  
    for (int j = n - 1; j >= iq + 1; --j) {
        /* The Givens rotation is done with the matrix (cc cs, cs -cc).
           If cc is one, then element (j) of d is zero compared with element
           (j - 1). Hence we don't have to do anything. 
           If cc is zero, then we just have to switch column (j) and column (j - 1) 
           of J. Since we only switch columns in J, we have to be careful how we
           update d depending on the sign of gs.
           Otherwise we have to apply the Givens rotation to these columns.
           The i - 1 element of d has to be updated to h. */
        cc = d[j - 1];
        ss = d[j];
        h = distance(cc, ss);
        if (std::abs(h) < EPS) // h == 0
            continue;
        d[j] = 0.0;
        ss = ss / h;
        cc = cc / h;
        if (cc < 0.0) {
            cc = -cc;
            ss = -ss;
            d[j - 1] = -h;
        }
        else
            d[j - 1] = h;
        xny = ss / (1.0 + cc);
        for (int k = 0; k < n; k++) {
            t1 = J(k, j - 1);
            t2 = J(k,j);
            J(k, j - 1) = t1 * cc + t2 * ss;
            J(k, j) = xny * (t1 + J(k, j - 1)) - t2;
        }
    }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
       into column iq - 1 of R
    */
    for (int i = 0; i < iq; i++)
        R(i, iq - 1) = d[i];
  
    if (std::abs(d[iq - 1]) <= EPS * R_norm) {
        // problem degenerate
        return false;
    }
    R_norm = std::max<double>(R_norm, std::abs(d[iq - 1]));
    return true;
}

void QPSolver::Impl::delete_constraint(int& iq, int l) {
    int qq = -1; // just to prevent warnings from smart compilers
    double cc, ss, h, xny, t1, t2;
  
    /* Find the index qq for active constraint l to be removed */
    for (int i = p; i < iq; ++i) {
        if (A[i] == l) {
            qq = i;
            break;
        }
    }
 
    /* remove the constraint from the active set and the duals */
    for (int i = qq; i < iq - 1; ++i) {
        A[i] = A[i + 1];
        u[i] = u[i + 1];
        R.col(i) = R.col(i + 1);
    }

    A[iq - 1] = A[iq];
    u[iq - 1] = u[iq];
    A[iq] = 0; 
    u[iq] = 0.0;
    R.row(iq - 1).setZero();

    /* constraint has been fully removed */
    iq--;
  
    if (iq == 0)
        return;
  
    for (int j = qq; j < iq; ++j) {
        cc = R(j, j);
        ss = R(j + 1, j);
        h = distance(cc, ss);
        if (std::abs(h) < EPS) // h == 0
            continue;
        cc = cc / h;
        ss = ss / h;
        R(j + 1, j) = 0.0;
        if (cc < 0.0) {
            R(j,j) = -h;
            cc = -cc;
            ss = -ss;
        } else {
            R(j,j) = h;
        }
    
        xny = ss / (1.0 + cc);
        for (int k = j + 1; k < iq; ++k) {
            t1 = R(j,k);
            t2 = R(j + 1, k);
            R(j,k) = t1 * cc + t2 * ss;
            R(j + 1, k) = xny * (t1 + R(j,k)) - t2;
        }
        for (int k = 0; k < n; ++k) {
            t1 = J(k,j);
            t2 = J(k, j + 1);
            J(k,j) = t1 * cc + t2 * ss;
            J(k, j + 1) = xny * (J(k,j) + t1) - t2;
        }
    }
}

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
    _data->g = g.asEigen();
    _data->c = c.asEigen();
    _data->isDiagonal = false;
    _data->maxIterations = 100;
}

QPSolver::QPSolver(
    lsst::ndarray::Array<double const,1,0> const & g,
    lsst::ndarray::Array<double const,1,0> const & c
) : _data(boost::make_shared<Data>())
{
    int n = c.getSize<0>();
    detail::checkSize(
        g.getSize<0>(), n, 
        "Size of QP matrix diagonal (%d) does not match vector size (%d)."
    );
    _data->g = g.asEigen().asDiagonal();
    _data->c = c.asEigen();
    _data->isDiagonal = true;
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
    _data->isDiagonal = false;
    _data->maxIterations = 100;
}

QPSolver::QPSolver(
    Eigen::VectorXd const & g,
    Eigen::VectorXd const & c
) : _data(boost::make_shared<Data>())
{
    int n = c.size();
    detail::checkSize(
        g.size(), n, 
        "Size of QP matrix diagonal (%d) does not match vector size (%d)."
    );
    _data->g = g.asDiagonal();
    _data->c = c;
    _data->isDiagonal = true;
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
    _data->ae = a.asEigen().transpose();
    _data->be = b.asEigen(); 
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
    _data->ai = a.asEigen().transpose();
    _data->bi = b.asEigen(); 
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
    double cost = Impl(*_data, qp_x).solve();
    detail::checkSize(
        x.getSize<0>(), qp_x.size(),
        "Size of solution vector (%d) is incorrect (%d)."
    );
    x.asEigen() = qp_x;
    return cost;
}

double QPSolver::solve(Eigen::VectorXd & x) const {
    return Impl(*_data, x).solve();
}

}}} // namespace lsst::meas::multifit
