// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018  AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "lsst/meas/modelfit/RegularizedMoments.h"

namespace lsst {
namespace meas {
namespace modelfit {

namespace {

using Moments = MomentsModel::Moments; 

std::pair<Moments, Moments> buildTestMoments(){
    Moments Q(6, 4, 3, 2, 1, 4);
    Moments W(2, 4, 3.1, 2.5, 1.2, 3.7);
    return std::make_tuple<Q, W>;
}

struct alphaX {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (Q.second.inverse() + W.second.inverse()).inverse()*
               (Q.second.inverse()*Q.first + W.second.inverse()*W.first)(0);
    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;
        vec(0, 0) = 0;
        vec(1, 0) = −1*(Q[5]*W[3]−Q[4]*W[4]−W[4]^2+W[3]*W[5])/
                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        vec(2, 0) = (Q[4]*W[3]−Q[3]*W[4])/
                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec3Top = W[2]*W[4]^3+W[1]*W[3]*W[5]^2−(Q[2]*W[4]−W[2]*W[4])*Q[4]^2−
                         (Q[1]*W[3]−W[1]*W[3])*Q[5]^2+(W[4]^2*W[5]-W[3]*W[5]^2)*Q[1]−
                         (W[4]^3−W[3]*W[4]*W[5])*Q[2]+(2*W[2]*W[4]^2+Q[1]*W[4]*W[5]−
                         (2*W[4]^2−W[3]*W[5])*Q[2]−(W[2]*W[3]+W[1]*W[4])*W[5])*Q[4]+
                         (Q[2]*W[3]*W[4]−W[2]*W[3]*W[4]−W[1]*W[4]^2+2*W[1]*W[3]*W[5]+
                         (W[4]^2−2*W[3]*W[5])*Q[1]+(Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−
                         W[1]*W[4])*Q[4])*Q[5]−(W[2]*W[3]*W[4]+W[1]*W[4]^2)*W[5];

        double vec3Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                            W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                            (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q3+
                            4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−2*((Q[3]+W[3])*Q[4]^2+
                            W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−2*W[3]*W[5])*Q[3]+
                            2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(3, 0) = vec3Top/vec3Bottom;

        double vec4Top = W[2]*W[3]*W[4]^2−W[1]*W[4]^3+(Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[4]^2+
                         (W[4]^3−W[3]*W[4]*W[5])*Q[1]−(W[3]*W[4]^2−W[3]^2*W[5])*Q[2]+
                         (2*W[2]*W[4]^2+Q[1]*W[4]*W[5]−(2*W[4]^2−W[3]*W[5])*Q[2]−
                         (W[2]*W[3]+W[1]*W[4])*W[5])*Q[3]−2*(W[1]*W[4]^2−W[1]*W[3]*W[5]−
                         (W[4]^2−W[3]*W[5])*Q[1]+(Q[2]*W[4]−W[2]*W[4])*Q[3])*Q[4]+
                         (Q[2]*W[3]^2−W[2]*W[3]^2−Q[1]*W[3]*W[4]+W[1]*W[3]*W[4]+
                         (Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[3]−2*(Q[1]*W[3]−
                         W[1]*W[3])*Q[4])*Q[5]−(W[2]*W[3]^2−W[1]*W[3]*W[4])*W[5];

        double vec4Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                            W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                            (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+
                            4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−2*((Q[3]+W[3])*Q[4]^2+
                            W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−2*W[3]*W[5])*Q[3]+
                            2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(4, 0) = -1*vec4Top/vec4Bottom;

        double vec5Top = (Q[2]*W[4]−W[2]*W[4])*Q[3]^2+(Q[1]*W[3]−W[1]*W[3])*Q[4]^2+
                         (Q[2]*W[3]*W[4]−W[2]*W[3]*W[4]−Q[1]*W[4]^2+W[1]*W[4]^2)*Q[3]−
                         (Q[2]*W[3]^2−W[2]*W[3]^2−Q[1]*W[3]*W[4]+W[1]*W[3]*W[4]+
                          (Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[3])*Q[4];

        double vec5Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+
                            Q[3]^2*W[5]^2+W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−
                            W[3]*W[5])*Q[4]^2+(Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−
                            2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+4*(W[4]^3−Q[3]*W[4]*W[5]−
                            W[3]*W[4]*W[5])*Q[4]−2*((Q[3]+W[3])*Q[4]^2+W[3]*W[4]^2−
                            Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−2*W[3]*W[5])*Q[3]+
                            2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(5, 0) = -1*vec5Top/vec5Bottom;
        return vec;
    }
};

struct alphaY {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (Q.second.inverse() + W.second.inverse()).inverse()*
               (Q.second.inverse()*Q.first + W.second.inverse()*W.first)(1);
    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;
        vec(0, 0) = 0;

        vec(1, 0) = (Q[4]*W[5] - Q[5]*W[4])/
                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        vec(2, 0) = (Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])/
                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec3Top = (Q[2]*W[5]−W[2]*W[5])*Q[4]^2+(Q[1]*W[4]−W[1]*W[4])*Q[5]^2+
                         (Q[2]*W[4]*W[5]−W[2]*W[4]*W[5]−Q[1]*W[5]^2+W[1]*W[5]^2)*Q[4]−
                         (Q[2]*W[4]^2−W[2]*W[4]^2−Q[1]*W[4]*W[5]+W[1]*W[4]*W[5]+
                         (Q[2]*W[4]−W[2]*W[4]+Q[1]*W[5]−W[1]*W[5])*Q[4])*Q[5];

        double vec3Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]W[4]^2*W[5]+
                            Q[3]^2*W[5]^2+W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                            (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+
                            4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−
                            2*((Q[3]+W[3])*Q[4]^2+W[3]W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+
                            (W[4]^2−2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(3, 0) = -1*vec3Top/vec3Bottom;

        vec4Top = W[2]*W[4]^3+W[1]*W[3]*W[5]^2−(Q[2]*W[4]−W[2]*W[4]+Q[1]*W[5]−W[1]*W[5])*Q[4]^2+
                  (W[4]^2*W[5]−W[3]*W[5]^2)*Q[1]−(W[4]^3−W[3]*W[4]*W[5])*Q[2]+
                  (Q[2]*W[4]*W[5]−W[2]*W[4]*W[5]−Q[1]*W[5]^2+W[1]*W[5]^2)*Q[3]+
                  2*(W[2]*W[4]^2−W[2]*W[3]*W[5]−(W[4]^2−W[3]*W[5])*Q[2]+
                  (Q[2]*W[5]−W[2]*W[5])*Q[3])*Q[4]−(Q[2]*W[3]*W[4]−W[2]*W[3]*W[4]+
                  2*W[1]W[4]^2−W[1]*W[3]*W[5]−(2*W[4]^2−W[3]*W[5])*Q1+(Q[2]*W[4]−
                  W[2]*W[4]+Q[1]*W[5]−W[1]*W[5])*Q[3]−2*(Q[1]*W[4]−W[1]*W[4])*Q[4])*Q[5]−
                  (W[2]*W[3]*W[4]+W[1]*W[4]^2)*W[5];

        vec4Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                     W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                     (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+
                     4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−2*((Q[3]+W[3])*Q[4]^2+
                     W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−2*W[3]*W[5])*Q[3]+
                     2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(4, 0) = vec4Top/vec4Bottom;

        vec5Top = W[2]*W[3]*W[4]^2−W[1]*W[4]^3+(Q[2]*W[5]−W[2]*W[5])*Q[3]^2+
                  (Q[1]*W[4]−W[1]*W[4])*Q[4]^2+(W[4]^3−W[3]*W[4]*W[5])*Q[1]−
                  (W[3]*W[4]^2−W[3]^2*W[5])*Q[2]+(W[2]*W[4]^2−Q[1]*W[4]*W[5]−
                  (W[4]^2−2*W[3]*W[5])*Q[2]−(2*W[2]*W[3]−W[1]*W[4])*W[5])*Q[3]−
                  (Q[2]*W[3]*W[4]−W[2]*W[3]*W[4]+2*W[1]*W[4]^2−W[1]*W[3]*W[5]−
                  (2*W[4]^2−W[3]*W[5])*Q[1]+(Q[2]*W[4]−W[2]*W[4]+Q[1]*W[5]−
                  W[1]*W[5])*Q[3])*Q[4]−(W[2]*W[3]^2−W[1]*W[3]*W[4])*W[5];

        vec5Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                     W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                     (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+
                     4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−2*((Q[3]+W[3])*Q[4]^2+
                     W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−2*W[3]*W[5])*Q[3]+
                     2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(5, 0) = -1*vec5Top/vec5Bottom;
        return vec;
    }
}

struct BetaX {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (Q[4]^2*W[3]−Q[3]*Q[5]W[3]+(W[4]^2−W[3]*W[5])*Q[3])/
               (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])
    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;

        vec(0, 0) = 0;
        vec(1, 0) = 0;
        vec(2, 0) = 0;

        double bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+W[3]^2*W[5]^2+
                        2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+(Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−
                        2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+4*(W[4]^3−Q[3]*W[4]W[5]−W[3]*W[4]W[5])*Q[4]−
                        2*((Q[3]+W[3])*Q[4]^2+W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−
                        2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        double vec3Top = Q[5]^2*W[3]^2+Q[4]^2*W[4]^2+W[4]^4−2*W[3]*W[4]^2*W[5]+W[3]^2*W[5]^2+
                         2*(W[4]^3−W[3]*W[4]W[5])*Q[4]−2*(Q[4]*W[3]W[4]+W[3]*W[4]^2−W[3]^2*W[5])*Q[5];

        vec(3, 0) = vec3Top/bottom;

        double vec4Top = 2*(Q[4]^2*W[3]*W[4]−(W[4]^3−W[3]*W[4]W[5])*Q[3]−(Q[3]*W[4]^2−W[3]*W[4]^2+
                         W[3]^2*W[5])*Q[4]−(Q[4]*W[3]^2−Q[3]*W[3]W[4])*Q[5]);

        vec(4, 0) = vec4Top/bottom;

        double vec5Top = Q[4]^2*W[3]^2−2*Q[3]*Q[4]W[3]*W[4]+Q[3]^2*W[4]^2;

        vec(5, 0) = vec5Top/bottom;

        return vec;
    }
};

struct BetaXY {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (Q[4]^2*W[4]−Q[3]*Q[5]W[4]+(W[4]^2−W[3]*W[5])*Q[4])/
               (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);
    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;

        vec(0, 0) = 0;
        vec(1, 0) = 0;
        vec(2, 0) = 0;

        double bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+W[3]^2*W[5]^2+
                        2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+(Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−
                        2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+4*(W[4]^3−Q[3]*W[4]W[5]−W[3]*W[4]W[5])*Q[4]−
                        2*((Q[3]+W[3])*Q[4]^2+W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−
                        2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        double vec3Top = Q[5]^2*W[3]*W[4]+Q[4]^2*W[4]*W[5]+(W[4]^2*W[5]−W[3]*W[5]^2)*Q[4]−
                         (W[4]^3−W[3]*W[4]W[5]+(W[4]^2+W[3]*W[5])*Q[4])*Q[5];

        vec(3, 0) = vec3Top/bottom;

        double vec4Top = W[4]^4−2*W[3]*W[4]^2*W[5]+W[3]^2*W[5]^2+(W[4]^2+W[3]*W[5])*Q[4]^2−(W[4]^2*W[5]−
                         W[3]*W[5]^2)*Q[3]+2*(W[4]^3−Q[3]*W[4]W[5]−W[3]*W[4]W[5])*Q[4]−(2*Q[4]*W[3]W[4]+
                         W[3]*W[4]^2−W[3]^2*W[5]−(W[4]^2+W[3]*W[5])*Q[3])*Q[5];

        vec(4, 0) = vec4Top/bottom;

        double vec5Top = Q[4]^2*W[3]*W[4]+Q[3]^2*W[4]*W[5]−(W[4]^3−W[3]*W[4]W[5])*Q[3]+(W[3]*W[4]^2−
                         W[3]^2*W[5]−(W[4]^2+W[3]*W[5])*Q[3])*Q[4];

        vec(5, 0) = vec5Top/bottom;

        return vec;
    }
};

struct BetaY {
    static double computeValue(Moments const & Q, Moments const & W) {
        return (Q[4]^2*W[5]+(W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[5])/
            (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])
    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;

        vec(0, 0) = 0;
        vec(1, 0) = 0;
        vec(2, 0) = 0;

        double bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+W[3]^2*W[5]^2+
                        2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+(Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−
                        2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+4*(W[4]^3−Q[3]*W[4]W[5]−W[3]*W[4]W[5])*Q[4]−
                        2*((Q[3]+W[3])*Q[4]^2+W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−
                        2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        double vec3Top = Q[5]^2*W[4]^2−2*Q[4]*Q[5]W[4]*W[5]+Q[4]^2*W[5]^2;

        vec(3, 0) = vec3Top/bottom;

        double vec4Top = 2*(Q[4]^2*W[4]*W[5]+(W[4]^2*W[5]−Q[3]*W[5]^2−W[3]*W[5]^2)*Q[4]−(Q[4]*W[4]^2+
                         W[4]^3−Q[3]*W[4]W[5]−W[3]*W[4]W[5])*Q[5]);

        vec(4, 0) = vec4Top/bottom;

        double vec5Top = Q[4]^2*W[4]^2+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+W[3]^2*W[5]^2−
                         2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+2*(W[4]^3−Q[3]*W[4]W[5]−W[3]*W[4]W[5])*Q[4];

        vec(5, 0) = vec5Top/bottom;

        return vec;
    }
};

struct Norm {
    static double computeValue(Moments Q, Moments W) {
        double expTop = Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+2*(W[2]*W[4]−
                        W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+(Q[2]^2−2*Q[2]*W[2]+
                        W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+
                        W[1]^2)*Q[5]
        
        double expBottom = 2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])

        double bottom = std::sqrt(−4*pi^2*(Q[4]+W[4])^2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]))

        return std::exp(expTop/expBottom)/bottom

    }

    static Moments computeGradient(Moments const & Q, Moments const & W) {
        Moments vec;

        vec(0,0) = 0;

        double vec1Top = ((Q[2]−W[2])*Q[4]−(Q[1]−W[1])*Q[5]+
                          Q[2]*W[4]−W[2]*W[4]−Q[1]*W[5]+W[1]*W[5])*
                         std::exp((Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                   2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                   (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                   W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])/
                                  (2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])));
        double vec1Bottom = std::sqrt(−4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]))(Q[4]^2−(Q[3]+W[3])*
                            Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        vec(1,0) = -1*vec1Top/vec1Bottom;

        double vec2Top = ((Q[2]−W[2])*Q[3]−(Q[1]−W[1])*Q[4]+
                          Q[2]*W[3]−W[2]*W[3]−Q[1]*W[4]+W[1]*W[4])*
                         std::exp((Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                   2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                   (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                   W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])/
                                  (2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])));

        double vec2Bottom = np.sqrt(−4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]))*
                            (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        vec(2,0) = vec2Top/vec2Bottom; 

        double vec3FirstTop = 2*pi^2*(Q[5]+W[5])*
                              std::exp((Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                     2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                     (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                     W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])/
                                     (2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])));

        double vec3FirstBottom = (−4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]))^(3/2);

        double vec3SecondTopFirst = (Q[2]^2−2*Q[2]*W[2]+W[2]^2)\
                                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec3SecondTopSecondTop = (Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                         2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                         (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                         W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])(Q[5]+W[5]);

        double vec3SecondTopSecondBottom = (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])^2;

        double vec3SecondExpTop = Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                  2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                  (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                  W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5];

        double vec3SecondExpBottom = 2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec3SecondBottom = 2*std::sqrt(−4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]));

        vec(3,0) = -1*vec3FirstTop/vec3FirstBottom +
                   ((vec3SecondTopFirst+(vec3SecondTopSecondTop)/(vec3SecondTopSecondBottom))*
                    std::exp(vec3SecondExpTop/vec3SecondExpBottom))/vec3SecondBottom;

        double vec4FirstTop = 4*pi^2*(Q[4]+W[4])*
                              std::exp((Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                        2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                        (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                        W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])/
                                       (2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])));
        
        double vec4FirstBottom = (−4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]))^(3/2);
        
        double vec4SecondTopFirst = ((Q[1]−W[1])*Q[2]−Q[1]*W[2]+W[1]*W[2])/
                                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec4SecondTopSecondTop = (Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                         2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                         (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                         W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])(Q[4]+W[4]);

        double vec4SecondTopSecondBottom = (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])^2;

        double vec4SecondTopSecondExpTop = std::exp(Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+
                                                    W[1]^2*W[5]+2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+
                                                    Q[1]*W[4]−W[1]*W[4])*Q[2]+(Q[2]^2−2*Q[2]*W[2]+
                                                    W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                                    W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5]);

        double vec4SecondTopSecondExpBottom = 2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec4SecondBottom = std::sqrt(-4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]));

        vec(4, 0) = vec4FirstTop/vec4FirstBottom - (vec4SecondTopFirst +
                                                    vec4SecondTopSecondTop/vec4SecondTopSecondBottom*
                                                    std::exp(vec4SecondTopSecondExpTop/
                                                             vec4SecondTopSecondExpBottom))/
                                                    vec4SecondBottom;

        double vec5FirstTop = 4*pi^2*(Q[4]+W[4])*
                              std::exp((Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                        2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                        (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                         W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5])/
                                       (2*(Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])));

        double vec5FirstBottom = (−4*pi^2*(Q[4]+W[4])2+4*pi^2*(Q[3]+W[3])(Q[5]+W[5]))^(3/2);

        double vec5SecondTopFirst = (Q[1]^2−2*Q[1]*W[1]+W[1]^2)/
                                    (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5]);

        double vec5SecondTopSecondTop = (Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]
                                         +2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−
                                         W[1]*W[4])*Q[2]+(Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−
                                         2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+W[1]*W[2])*Q[4]+(Q[1]^2−
                                         2*Q[1]*W[1]+W[1]^2)*Q[5])(Q[3]+W[3]);

        double vec5SecondTopSecondBottom = (Q[4]^2−(Q[3]+W[3])*Q[5]+2*Q[4]*W[4]+W[4]^2−Q[3]*W[5]−W[3]*W[5])^2;

        double vec5SecondTopSecondExpTop = Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+W[1]^2*W[5]+
                                           2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[2]+
                                           (Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+
                                           W[1]*W[2])*Q[4]+(Q[1]^2−2*Q[1]*W[1]+W[1]^2)*Q[5];

        double vec5SecondTopSecondExpBottom = Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+
                                              W[1]^2*W[5]+2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+
                                              Q[1]*W[4]−W[1]*W[4])*Q[2]+(Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−
                                              2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+W[1]*W[2])*Q[4]+(Q[1]^2−
                                              2*Q[1]*W[1]+W[1]^2)*Q[5];

        double vec5SecondBottom = 2*std::sqrt(Q[2]^2*W[3]+W[2]^2*W[3]−2*W[1]*W[2]W[4]+Q[1]^2*W[5]+
                                              W[1]^2*W[5]+2*(W[2]*W[4]−W[1]*W[5])*Q[1]−2*(W[2]*W[3]+
                                              Q[1]*W[4]−W[1]*W[4])*Q[2]+(Q[2]^2−2*Q[2]*W[2]+W[2]^2)*Q[3]−
                                              2*((Q[1]−W[1])*Q[2]−Q[1]*W[2]+W[1]*W[2])*Q[4]+(Q[1]^2−
                                              2*Q[1]*W[1]+W[1]^2)*Q[5]);

        vec(5, 0) = -1*vec5FistTop/vec5FistBottom +
                    (vec5SecondTopFirst+(vec5SecondTopSecondTop/vec5SecondTopSecondBottom)*
                    std::exp(vec5SecondTopSecondExpTop/vec5SecondTopSecondExpBottom))/
                    vec5SecondBottom;

        return vec;

    }
private:
    using afw::geom::pi

};

MomentsModel::FirstMoment makeAlpha(Moments const & Q, Moments const & W) {
    double x = AlphaX::computeValue(Q, W);
    double y = AlphaY::computeValue(Q, W);

    return MomentsModel::FirstMoment(x, y);
};

auto makeAlphaGrad(Moments const & Q, Moments const & W) {
    x = AlphaX::computeGradient(Q, W);
    y = AlphaY::computeGradient(Q, W);

    return make_tuple(x, y);
}



MomentsModel::SecondMoment makeBeta(Moments const & Q, Moments const & W) {
    double x = BetaX::computeValue(Q, W);
    double xy = BetaXY::computeValue(Q, W);
    double y = BetaY::computeValue(Q, W);

    MomentsModel::SecondMoment beta;
    
    beta(0, 0) = x;
    beta(1, 0) = xy;
    beta(0, 1) = xy;
    beta(1, 1) = y;

    return y
}

auto makeBetaGrad(Moments const & Q, Moments const & W) {
    return std::make_tuple(BetaX::computeGradient(Q, W), BetaXY::computeGradient(Q, W),
                           BetaY::computeGadient(Q, W));
}

bool approxEqual(Moments const & first, Moments const & second, double tol=1e-6){
    for (size_t i = 0; i < 6; ++i){
        if (abs(first(i, 0) - second(i, 0)) > tol) {
            return false;
        }
    }
    return true;
}

} //end anonymous

bool testAlphaX(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = AlphaX.computeValue(Q, W);
    Moments firstRes = AlphaX.computeGradient(Q, W);
    zeroTruth = 4.00033545790003;
    Moments firstTruth({0,
                        0.557195571955720,
                        −0.00335457900033546,
                        −0.00411214444247764,
                        0.00843596158202444,
                        −0.0000506394012127103});
    if (abs(zeroTruth - zeroRes) > tol) {
        return false;
    }
    return approxEqual(firstMoments, firstTruth, tol);
}

bool testAlphaY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = AlphaY.computeValue(Q, W);
    Moments firstRes = AlphaY.computeGradient(Q, W);
    double zeroTruth = 3.05300234820530;
    Moments firstTruth({0,
                        0.0369003690036900,
                        0.469976517946998,
                        −0.000272327446521693,
                        −0.00291142797372289,
                        0.00709458010990103});

    if (abs(zeroTruth - zeroRes) > tol) {
        return false
    }
    return approxEqual(firstMoments, firstTruth, tol)
}

bool testBetaX(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaX.computeValue(Q, W);
    Moments firstRes = BetaX.computeGradient(Q, W);
    double zeroTruth = 1.11103656491110;
    Moments firstTruth({0,
                        0,
                        0,
                        0.310466905407061,
                        −0.00373831312952513,
                        0.0000112532002694914});

    if (abs(zeroTruth - zeroRes) > tol) {
        return false
    }
    return approxEqual(firstMoments, firstTruth, tol)
}

bool testBetaXY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaXY.computeValue(Q, W);
    Moments firstRes = BetaXY.computeGradient(Q, W);
    double zeroTruth = 0.543777255954378;
    Moments firstTruth({0,
                        0,
                        0,
                        0.0205607222123882,
                        0.261745049520270,
                        −0.00157657335775577});

    if (abs(zeroTruth - zeroRes) > tol) {
        return false
    }
    return approxEqual(firstMoments, firstTruth, tol)
}

bool testBetaY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaY.computeValue(Q, W);
    Moments firstRes = BetaY.computeGradient(Q, W);
    double zeroTruth = 1.91680644079168;
    Moments firstTruth({0,
                        0,
                        0,
                        0.00136163723260849,
                        0.0346846138706271,
                        0.220877927421585});

    if (abs(zeroTruth - zeroRes) > tol) {
        return false
    }
    return approxEqual(firstMoments, firstTruth, tol)
}


bool testNorm(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = norm.computeValue(Q, W);
    Moments firstRes = norm.computeGradient(Q, W);
    double zeroTruth = 0.0915084542604366/afw::geom::pi;
    Moments firstTruth({0,
                        -0.000675339145833491/afw::geom::pi,
                        0.00138137552556848/afw::geom::pi,
                        −0.0118159430257175/afw::geom::pi,
                        0.00674319680500958/afw::geom::pi,
                        0.00689645127785071/afw::geom::pi});

    if (abs(zeroTruth - zeroRes) > tol) {
        return false
    }
    return approxEqual(firstMoments, firstTruth, tol)
}

void MomentsModel::at(Moments const & inputQ) {
   Q = Moments(inputQ); 
   makeValue();
}

void MomentsModel::makeValue() {
    alpha = makeAlpha(Q, W);
    beta = makeBeta(Q, W);
    norm = Norm::computeValue(Q, W);

    double zero = 2*afw::geom::pi*Q(0, 0)*beta.determinant()*norm;
    FistMoment one = zero*alpha;
    SecondMoment two = zero*(beta + alpha*alpha.transpose());

    value = Moments(zero, one(0, 0), one(1, 0), two(0, 0), two(0, 1), two(1, 1));
}

Moments MomentsModel::computeValue() {
    return value;
}

Jacobian MomentsModel::computeJacobian() {
    // Calculate the matrices that will be needed for the rest of this
    // function
    Moments normGrad = Norm::computeGradient(Q, W);

    Moments alphaXGrad, alphaYGrad;
    std::tie(alphaXGrad, alphaYGrad;

    Moments betaXGrad, betaXYGrad, betaYGrad;
    std::tie(betaXGrad, betaXYGrad, betaYGrad) = makeBetaGrad(Q, W);

    // Calculate the gradient along the first moment
    Moments zerothGrad;
    for (size_t i=0; i < 5; ++i) {
        // the first term in the equation is only needed for the derivative
        // with respect to Q_0, otherwise the term always is zero
        float accumulator = 0;
        if (i == 0) {
            accumulator += beta.determinant()*norm;
        }

        accumulator += Q[0]*norm*(betaXGrad[i]*beta(1,1) + beta(0,0)*betaYGrad[i] - 2*beta(1,0)*betaXY[i]);

        accumulator += Q[0]*beta.determinant()*normGrad[i];

        zerothGrad(i,0) = 2*afw::geom::pi*accumulator;
    }

    // Calculate the gradient along the two components of the first moment
    Moments firstX, firstY;

    firstX = alpha(0, 0)*zerothGrad + value(0, 0)*alphaXGrad;
    firstY = alpha(1, 0)*zerothGrad + value(0, 0)*alphaYGrad;

    // Calculate the gradient along each of the components of the second
    // moment
    SecondMoment modBeta = beta + alpha*alpha.transpose();

    Moments secondX, secondXY, secondY;

    secondX = modBeta(0, 0)*normGrad + value(0, 0)*(betaXGrad + 2*alpha(0, 0)*alphaXGrad);
    secondXY = modBeta(0, 1)*normGrad + 
               value(0, 0)*(betaXYGrad + 2*(alphaXGrad*alpha(1, 0) +alphaYGrad*alpha(0, 0)));
    secondX = modBeta(1, 1)*normGrad + value(0, 0)*(betaYGrad + 2*alpha(1, 0)*alphaYGrad);

    // Build the result and return it
    Jacobian result;

    result(0) = zerothGrad.transpose();
    result(1) = firstX.transpose();
    result(2) = firstY.transpose();
    result(3) = secondX.transpose();
    result(4) = secondXY.transpose();
    result(5) = secondY.transpose();

    return result;
}
}}} // Close lsst::meas::modelfit
