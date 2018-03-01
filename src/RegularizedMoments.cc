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

struct alphaX {
    static double computeValue(Moments Q, Moments W) {
        return (Q.second.inverse() + W.second.inverse()).inverse()*
               (Q.second.inverse()*Q.first + W.second.inverse()*W.first)(0)
    }

    static Moments::ParameterVector computeGradient(Moments Q, Moments W) {
        Moments::ParameterVector vec;
        vec(0, 0) = 0;
        vec(1, 0) = W[4]^2−W[3]*W[5]+W[4]*Q[4]−W[3]Q[4]/
                    (W[4]^2−W[3]*W[5]−W[5]*Q[3]+2*W[4]*Q[4]+Q[4]^2−(W[3]+Q[3])*Q[5]);
        vec(2, 0) = −W[4]*Q[3]−W[3]*Q[4]/(W[4]^2−W[3]*W[4]−W[4]*Q[3]+2*W[4]*Q[4]+Q[4]^2−(W[3]+Q[3])*Q[5]);

        double vec3Top = W[2]*W[4]^3+W[1]*W[3]*W[5]^2−(Q[2]*W[4]−W[2]*W[4])*Q[4]^2−
                         (Q[1]*W[3]−W[1]*W[3])*Q[5]^2+(W[4]^2*W[5]−W[3]*W[5]^2)*Q[1]−
                         (W[4]^3−W[3]*W[4]*W[5])*Q[2]+(2*W[2]*W[4]^2+Q[1]*W[4]*W[5]−
                         (2*W[4]^2−W[3]*W[5])*Q[2]−(W[2]*W[3]+W[1]*W[4])*W[5])*Q[4]+
                         (Q[2]*W[3]*W[4]−W[2]*W[3]*W[4]−W[1]*W[4]^2+2*W[1]*W[3]*W[5]+
                         (W[4]^2−2*W[3]*W[5])*Q[1]+(Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−
                         W[1]*W[4])*Q[4])*Q[5]−(W[2]*W[3]*W[4]+W[1]*W[4]^2)*W[5];

        double vec3Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                            W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                            (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−W[3]*W[5]^2)*Q[3]+
                            4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−2*((Q[3]+W[3])*Q[4]^2+
                            W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+(W[4]^2−2*W[3]*W[5])*Q[3]+
                            2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(3, 0) = vec3Top/vec3Bottom;

        double vec4Top = W[2]*W[3]*W[4]^2−W[1]*W[4]^3+(Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−
                         W[1]*W[4])*Q[4]^2+(W[4]^3−W[3]*W[4]*W[5])*Q[1]−
                         (W[3]*W[4]^2−W[3]^2*W[5])*Q[2]+(2*W[2]*W[4]^2+Q[1]*W[4]*W[5]−
                         (2*W[4]^2−W[3]*W[5])*Q[2]−(W[2]*W[3]+W[1]*W[4])*W[5])*Q[3]−
                         2*(W[1]*W[4]^2−W[1]*W[3]*W[5]−(W[4]^2−W[3]*W[5])*Q[1]+
                         (Q[2]*W[4]−W[2]*W[4])*Q[3])*Q[4]+(Q[2]*W[3]^2−W[2]*W[3]^2−
                         Q[1]*W[3]*W[4]+W[1]*W[3]*W[4]+(Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−
                         W[1]*W[4])*Q[3]−2*(Q[1]*W[3]−W[1]*W[3])*Q[4])*Q[5]−(W[2]*W[3]^2−
                         W[1]*W[3]*W[4])*W[5];

        double vec4Bottom = Q[4]^4+4*Q[4]^3*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                            W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                            (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−
                            W[3]*W[5]^2)*Q[3]+4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−
                            2*((Q[3]+W[3])*Q[4]^2+W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+
                            (W[4]^2−2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(4, 0) = -1*vec4Top/vec4Bottom;

        double vec5Top = (Q[2]*W[4]−W[2]*W[4])*Q[3]^2+(Q[1]*W[3]−W[1]*W[3])*Q[4]^2+
                         (Q[2]*W[3]*W[4]−W[2]*W[3]*W[4]−Q[1]*W[4]^2+W[1]*W[4]^2)*Q[3]−
                         (Q[2]*W[3]^3−W[2]*W[3]^2−Q[1]*W[3]*W[4]+W[1]*W[3]*W[4]+
                         (Q[2]*W[3]−W[2]*W[3]+Q[1]*W[4]−W[1]*W[4])*Q[3])*Q[4];

        double vec5Bottom = Q[4]^4+4*Q[3]*4*W[4]+W[4]^4−2*W[3]*W[4]^2*W[5]+Q[3]^2*W[5]^2+
                            W[3]^2*W[5]^2+2*(3*W[4]^2−Q[3]*W[5]−W[3]*W[5])*Q[4]^2+
                            (Q[3]^2+2*Q[3]*W[3]+W[3]^2)*Q[5]^2−2*(W[4]^2*W[5]−
                            W[3]*W[5]^2)*Q[3]+4*(W[4]^3−Q[3]*W[4]*W[5]−W[3]*W[4]*W[5])*Q[4]−
                            2*((Q[3]+W[3])*Q[4]^2+W[3]*W[4]^2−Q[3]^2*W[5]−W[3]^2*W[5]+
                            (W[4]^2−2*W[3]*W[5])*Q[3]+2*(Q[3]*W[4]+W[3]*W[4])*Q[4])*Q[5];

        vec(5, 0) = -1*vec5Top/vec5Bottom;
        return vec;
    }
};

std::pair<Moments, Moments> buildTestMoments(){
    Moments Q(6, 4, 3, 2, 4, 1);
    Moments W(2, 4, 3.1, 2.5, 3.7, 1.2);
    return std::make_pair<Q, W>;
}

bool testAlphaX() {
    Moments Q, W;
    bool result = true;
    std::tie(Q, W) = buildTestMoments();
    Moments::ParameterVector zeroRes = AlphaX.computeValue(Q, W)
    Moments::ParameterVector firstRes = AlphaX.computeGradient(Q, W)
    zeroTruth =
}


struct alphaY {
    static double computeValue(Moments Q, Moments W) {
        return (Q.second.inverse() + W.second.inverse()).inverse()*
               (Q.second.inverse()*Q.first + W.second.inverse()*W.first)(1)
    }

    static Moments::ParameterVector computeGradient(Moments Q, Moments W) {
        Moments::ParameterVector vec;

    }
}

struct beta {
    static double computeValue(Moments Q, Moments W) {

    }

    static Moments::ParameterVector computeGradient(Moments Q, Moments W) {

    }
};

struct Norm {
    static double computeValue(Moments Q, Moments W) {

    }

    static Moments::ParameterVector computeGradient(Moments Q, Moments W) {

    }
};
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
    zero = *begin;
    // The next two entries will be the first moment, i.e. x, y
    // vec[1], vec[2]
    checkIter(begin++);
    double firstX = *begin;
    checkIter(begin++);
    double firstY = *begin;
    one << firstX, firstY;
    // The next three entries correspond to the second moment
    // i.e. xx, xy, yy: vec[3], vec[4], vec[5]
    checkIter(begin++);
    double secondX = *begin;
    checkIter(begin++);
    double secondXY = *begin;
    checkIter(begin++);
    double secondY = *begin;
    two << secondX, secondXY, secondY
}

Moments::Moments(std::vector<double> moments) {
    Moments(moments.begin(), moments.end());
}

Moments::Moments(std::initializer_list<double> initList) {
    Moments(initList.begin(), initList.end());
}

Moments::Moments(double zero, FirstMoment first, SecondMoment second): zero(zero), one(first), two(second){
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

}}} // Close lsst::meas::modelfit
