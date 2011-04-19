// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_RandomEngine
#define LSST_MEAS_MULTIFIT_RandomEngine

#include "lsst/ndarray.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

namespace lsst { namespace meas { namespace multifit {

class RandomEngine {
public:

    double drawNormal() { return _normal(); }

    double drawUniform() { return _uniform(); }

    RandomEngine();

    RandomEngine(RandomEngine const & other);

    RandomEngine & operator=(RandomEngine const & other);

private:
    typedef boost::mt19937 Engine;
    typedef boost::variate_generator< Engine &, boost::normal_distribution<double> > NormalGenerator;
    typedef boost::variate_generator< Engine &, boost::uniform_real<double> > UniformGenerator;
    
    Engine _engine;
    NormalGenerator _normal;
    UniformGenerator _uniform;
};



}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_RandomEngine
