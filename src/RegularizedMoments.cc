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

    static Moments::ParameterVector computeGradient(Moments const & Q, Moments const & W) {
        Moments::ParameterVector vec;
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

    static Moments:ParameterVector computeGradient(Moments const & Q, Moments const & W) {
        Moments::ParameterVector vec;
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

    static Moments::ParameterVector computeGradient(Moments const & Q, Moments const & W) {
        Moments::ParameterVector vec;

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

    static Moments::ParameterVector computeGradient(Moments const & Q, Moments const & W) {
        Moments::ParameterVector vec;

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

    static Moments::ParameterVector computeGradient(Moments const & Q, Moments const & W) {
        Moments::ParameterVector vec;

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

    static Moments::ParameterVector computeGradient(Moments const & Q, Moments const & W) {
        Moments::ParameterVector vec;

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

Moments makeBeta(Moments const & Q, Moments const & W) {
    double x = BetaX::computeValue(Q, W);
    double xy = BetaXY::computeValue(Q, W);
    double y = BetaY::computeValue(Q, W);

    Moments::SecondMoment beta;
    
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

} //end anonymous

bool testAlphaX(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = AlphaX.computeValue(Q, W);
    Moments::ParameterVector firstRes = AlphaX.computeGradient(Q, W);
    Moments resMoments(firstRes);
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
    return firstMoments.aproxEqual(firstTruth, tol);
}

bool testAlphaY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = AlphaY.computeValue(Q, W);
    Moments::ParameterVector firstRes = AlphaY.computeGradient(Q, W);
    Moments resMoments(firstRes);
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
    return firstMoments.aproxEqual(firstTruth, tol)
}

bool testBetaX(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaX.computeValue(Q, W);
    Moments::ParameterVector firstRes = BetaX.computeGradient(Q, W);
    Moments resMoments(firstRes);
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
    return firstMoments.aproxEqual(firstTruth, tol)
}

bool testBetaXY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaXY.computeValue(Q, W);
    Moments::ParameterVector firstRes = BetaXY.computeGradient(Q, W);
    Moments resMoments(firstRes);
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
    return firstMoments.aproxEqual(firstTruth, tol)
}

bool testBetaY(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = BetaY.computeValue(Q, W);
    Moments::ParameterVector firstRes = BetaY.computeGradient(Q, W);
    Moments resMoments(firstRes);
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
    return firstMoments.aproxEqual(firstTruth, tol)
}


bool testNorm(double tol) {
    Moments Q, W;
    std::tie(Q, W) = buildTestMoments();
    double zeroRes = norm.computeValue(Q, W);
    Moments::ParameterVector firstRes = norm.computeGradient(Q, W);
    Moments resMoments(firstRes);
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
    return firstMoments.aproxEqual(firstTruth, tol)
}

static double ZerothShapeMoment::computeMoment(Moments const & Q, Moments const & W){
    Moments::SecondMoment beta = makeBeta(Q, W);

    double norm = Norm::computeValue(Q, W);

    return 2*afw::geom::pi*Q[0]*beta.determinant()*norm;
}

static Moments::ParameterVector ZerothShapeMoment::computeGradient(Moments const & Q, Moments const & W){
    double norm = Norm::computeValue(Q, W);
    Moments::ParameterVector normGrad = Norm::computeGradient(Q, W);

    Moments::SecondMoment beta = makeBeta(Q, W);
    Moments::ParameterVector betaXGrad, betaXYGrad, betaYGrad;
    std::tie(betaXGrad, betaXYGrad, betaYGrad) = makeBetaGrad(Q, W);

    Moments::ParameterVector vec;
    for (size_t i=0; i < 5; ++i) {
        // the first term in the equation is only needed for the derivative
        // with respect to Q_0, otherwise the term always is zero
        float accumulator = 0;
        if (i == 0) {
            accumulator += beta.determinant()*norm;
        }

        accumulator += Q[0]*(betaXGrad[i]*beta(1,1) + beta(0,0)*betaYGrad[i] - 2*beta(1,0)*betaXY[i]);

        accumulator += Q[0]*beta.determinant()*normGrad[i];

        vec(i,0) = 2*afw::geom::pi*accumulator;
    }
    return vec;
}


template <typename iter>
Moments::Moments(iter begin, iter end) {
    auto checkIter = [&end](iter & begin) {
        std::string errMessage("Input vector too short")
        if (begin == end) {
            throw std::length_error(errMessage);
        }
    }
    // The first entry in the vector will be zeroth moment
    // vec[0]
    zeroth = *begin;
    // The next two entries will be the first moment, i.e. x, y
    // vec[1], vec[2]
    checkIter(begin++);
    double firstX = *begin;
    checkIter(begin++);
    double firstY = *begin;
    first << firstX, firstY;
    // The next three entries correspond to the second moment
    // i.e. xx, xy, yy: vec[3], vec[4], vec[5]
    checkIter(begin++);
    double secondX = *begin;
    checkIter(begin++);
    double secondXY = *begin;
    checkIter(begin++);
    double secondY = *begin;
    second << secondX, secondXY, secondY
}

Moments::Moments(std::vector<double> moments) {
    Moments(moments.begin(), moments.end());
}

Moments::Moments(std::initializer_list<double> initList) {
    Moments(initList.begin(), initList.end());
}

Moments::Moments(double zero, FirstMoment first, SecondMoment second): zero(zero), one(first), two(second){
}

Moments::Moments(ParameterVector params) {
    zeroth = params(0, 0);

    first << params(1, 0), params(2, 0);

    second << params(3, 0), params(4, 0), params(4, 0), params(5, 0);
}


double Moments::operator[](int i){
    switch(i) {
        case 0: return zeroth;
        case 1: return first(0, 0);
        case 2: return first(1, 0);
        case 3: return second(0, 0);
        case 4: return second(0, 1);
        case 5: return second(1, 1);
    }
}

Moments::ParameterVector Moments::getParameterVector(){
    Moments::ParameterVector vec;
    vec(0, 0) = zeroth;
    vec(1, 0) = first(1, 0);
    vec(2, 0) = first(2, 0);
    vec(3, 0) = second(0, 0);
    vec(4, 0) = second(0, 1);
    vec(5, 0) = second(1, 1);
    return vec
}

bool Moments:aproxEqual(Moments const & other, double tol) {
    for (int i = 0; i < 6; ++i){
        if (abs(*this[i] - other[i]) > tol) {
            return false
        }
    }
    return true
}


}}} // Close lsst::meas::modelfit
