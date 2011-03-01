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

#ifndef LSST_MEAS_MULTIFIT_SAMPLING_IterativeImportanceSampler
#define LSST_MEAS_MULTIFIT_SAMPLING_IterativeImportanceSampler

#include "lsst/meas/multifit/sampling/MixtureDistribution.h"

#include <vector>

namespace lsst { namespace meas { namespace multifit {

class IterativeImportanceSampler {
public:

    int getIterationCount() const { return _samples.size(); }

    Table const & getTable(int n) const { return _samples[n]; }

    BaseEvaluator::Ptr getEvaluator() const;

    void run(int size);

private:
    BaseEvaluator::Ptr _evaluator;
    MixtureDistribution _proposal;
    std::vector<Table> _samples;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SAMPLING_IterativeImportanceSampler
