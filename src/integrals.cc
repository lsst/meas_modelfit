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

#include "boost/math/special_functions/erf.hpp"

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/integrals.h"

namespace lsst { namespace meas { namespace multifit { namespace detail {

// translation of matlab 'bvn.m' routines by Alan Genz:
// http://www.math.wsu.edu/faculty/genz/homepage

double phid(double z) {
    return 0.5*boost::math::erfc(-z / M_SQRT2);
}

double bvnu(double h, double k, double rho) {
    pex::logging::Debug log("meas.multifit.integrals");
    log.debug<8>("Starting bvnu: h=%g, k=%g, rho=%g", h, k, rho);
    if (h == std::numeric_limits<double>::infinity() || h == std::numeric_limits<double>::infinity()) {
        return 0.0;
    } else if (h == -std::numeric_limits<double>::infinity()) {
        if (k == -std::numeric_limits<double>::infinity()) {
            return 1.0;
        } else {
            return phid(-k);
        }
    } else if (k == -std::numeric_limits<double>::infinity()) {
        return phid(-h);
    } else if (rho == 0.0) {
        return phid(-h) * phid(-k);
    }
    Eigen::ArrayXd w0;
    Eigen::ArrayXd x0;
    Eigen::ArrayXd w;
    Eigen::ArrayXd x;
    // setup Gauss-Legendre quadrature points and weights
    if (std::abs(rho) < 0.3) {
        w0.resize(3);
        x0.resize(3);
        w.resize(6);
        x.resize(6);
        w0 << 0.1713244923791705, 0.3607615730481384, 0.4679139345726904;
        x0 << 0.9324695142031522, 0.6612093864662647, 0.2386191860831970;
    } else if (std::abs(rho) < 0.75) {
        w0.resize(6);
        x0.resize(6);
        w.resize(12);
        x.resize(12);
        w0 << 0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
            0.2031674267230659, 0.2334925365383547, 0.2491470458134029;
        x0 << 0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
            0.5873179542866171, 0.3678314989981802, 0.1252334085114692;
    } else {
        w0.resize(10);
        x0.resize(10);
        w.resize(20);
        x.resize(20);
        w0 << .01761400713915212, 0.04060142980038694, 0.06267204833410906,
            .08327674157670475, 0.1019301198172404, 0.1181945319615184,
            0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
            0.1527533871307259;
        x0 << 0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
            0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
            0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
            0.07652652113349733;
    }
    w << w0, w0;
    x << 1.0 - x0, 1.0 + x0;
    double hk = h*k;
    double bvn = 0.0;
    if (std::abs(rho) < 0.925) {
        double hs = 0.5*(h*h + k*k);
        double asr = 0.5*std::asin(rho);
        Eigen::ArrayXd sn = (asr*x).sin();
        bvn = w.matrix().dot(((sn*hk - hs)/(1.0 - sn.square())).exp().matrix());
        bvn = 0.5*bvn*asr/M_PI + phid(-h)*phid(-k);
    } else {
        if (rho < 0) {
            k = -k;
            hk = -hk;
        }
        if (std::abs(rho) < 1) {
            double as = 1 - rho*rho;
            double a = std::sqrt(as);
            double bs = (h - k)*(h - k);
            double asr = -0.5*(bs/as + hk);
            double c = (4.0 - hk)/8.0;
            double d = (12.0 - hk)/80.0;
            if (asr > -100.0) {
                bvn = a*std::exp(asr)*(1.0 - c*(bs - as)*(1.0 - d*bs)/3.0 + c*d*as*as);
            }
            if (hk > -100.0) {
                double b = std::sqrt(bs);
                double sp = std::sqrt(2.0*M_PI)*phid(-b/a);
                bvn = bvn - std::exp(-0.5*hk)*sp*b*(1.0 - c*bs*(1.0 - d*bs)/3.0);
            }
            a = 0.5*a;
            Eigen::ArrayXd xs1 = (a*x).square();
            Eigen::ArrayXd asr1 = -(bs/xs1 + hk)/2;
            Eigen::Array<bool,Eigen::Dynamic,1> xi = asr1 > -100;
            Eigen::ArrayXd asr2 = Eigen::ArrayXd::Zero(xi.cast<int>().sum());
            Eigen::ArrayXd xs2 = Eigen::ArrayXd::Zero(asr2.size());
            Eigen::ArrayXd w2 = Eigen::ArrayXd::Zero(asr2.size());
            for (int i = 0, j = 0; i < x.size(); ++i) {
                if (xi[i]) {
                    asr2[j] = asr1[i];
                    xs2[j] = xs1[i];
                    w2[j] = w[i];
                    ++j;
                }
            }
            Eigen::ArrayXd sp = 1.0 + c*xs2*(1.0 + 5.0*d*xs2);
            Eigen::ArrayXd rs = (1.0 - xs2).sqrt();
            Eigen::ArrayXd ep = (-0.5*hk*xs2/(1.0 + rs).square()).exp()/rs;
            bvn = (a*(w2.matrix().dot((asr2.exp()*(sp - ep)).matrix())) - bvn)/(2.0*M_PI);
        }
        if (rho > 0.0) {
            bvn = bvn + phid(-std::max(h, k));
        } else if (h >= k) {
            bvn = -bvn;
        } else {
            double l = (h < 0) ? (phid(k) - phid(h)) : (phid(-h) - phid(-k));
            bvn = l - bvn;
        }
    }
    return std::max(0.0, std::min(1.0, bvn));
}

}}}} // namespace lsst::meas::multifit::detail
