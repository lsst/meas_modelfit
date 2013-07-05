// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
 */

#include <Eigen/Eigenvalues>

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"

#include "lsst/meas/multifit/optimizer.h"

namespace lsst { namespace meas { namespace multifit {

Eigen::VectorXd solveTrustRegion(
    Eigen::MatrixXd const & F, Eigen::VectorXd const & g, double r, double tolerance
) {
    static double const ROOT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    static int const ITER_MAX = 10;
    pex::logging::Debug log("meas.multifit.optimizer.solveTrustRegion");
    double const r2 = r*r;
    double const r2min = r2 * (1.0 - tolerance) * (1.0 - tolerance);
    double const r2max = r2 * (1.0 + tolerance) * (1.0 + tolerance);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(F);
    double const threshold = ROOT_EPS * eigh.eigenvalues()[g.size() - 1];
    Eigen::VectorXd qtg = eigh.eigenvectors().adjoint() * g;
    Eigen::VectorXd tmp(g.size());
    Eigen::VectorXd x(g.size());
    double mu = 0.0;
    double xsn = 0.0;
    if (eigh.eigenvalues()[0] >= threshold) {
        log.debug<7>("Starting with full-rank matrix");
        tmp = (eigh.eigenvalues().array().inverse() * qtg.array()).matrix();
        x = -eigh.eigenvectors() * tmp;
        xsn = x.squaredNorm();
        if (xsn <= r2max) {
            log.debug<7>("Ending with unconstrained solution");
            // unconstrained solution is within the constraint; no more work to do
            return x;
        }
    } else {
        mu = -eigh.eigenvalues()[0] + 2.0*ROOT_EPS*eigh.eigenvalues()[g.size() - 1];
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        int n = 0;
        while (eigh.eigenvalues()[++n] < threshold);
        log.debug<7>("Starting with %d zero eigenvalue(s)", n);
        if ((qtg.head(n).array() < ROOT_EPS * g.lpNorm<Eigen::Infinity>()).all()) {
            x = -eigh.eigenvectors().rightCols(n) * tmp.tail(n);
            xsn = x.squaredNorm();
            if (xsn < r2min) {
                // Nocedal and Wright's "Hard Case", which is actually
                // easier: Q_1^T g is zero (where the columns of Q_1
                // are the eigenvectors that correspond to the
                // smallest eigenvalue \lambda_1), so \mu = -\lambda_1
                // and we can add a multiple of any column of Q_1 to x
                // to get ||x|| == r.  If ||x|| > r, we can find the
                // solution with the usual iteration by increasing \mu.
                double tau = std::sqrt(r*r - x.squaredNorm());
                x += tau * eigh.eigenvectors().col(0);
                log.debug<7>("Ending; Q_1^T g == 0, and ||x|| < r");
                return x;
            }
            log.debug<7>("Continuing; Q_1^T g == 0, but ||x|| > r");
        } else {
            x = -eigh.eigenvectors() * tmp;
            xsn = x.squaredNorm();
            log.debug<7>("Continuing; Q_1^T g != 0");
        }
    }
    int nIter = 0;
    while ((xsn < r2min || xsn > r2max) && ++nIter < ITER_MAX) {
        log.debug<7>("Iterating at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
        mu += xsn*(std::sqrt(xsn) / r - 1.0)
            / (qtg.array().square() / (eigh.eigenvalues().array() + mu).cube()).sum();
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        x = -eigh.eigenvectors() * tmp;
        xsn = x.squaredNorm();
    }
    log.debug<7>("Ending at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
    return x;
}

}}} // namespace lsst::meas::multifit
