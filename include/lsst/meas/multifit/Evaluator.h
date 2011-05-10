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

#ifndef LSST_MEAS_MULTIFIT_Evaluator
#define LSST_MEAS_MULTIFIT_Evaluator

#include "lsst/base.h"
#include "lsst/meas/multifit/definition/Definition.h"
#include "lsst/meas/multifit/grid/Grid.h"
#include "lsst/meas/multifit/BaseEvaluator.h"

namespace lsst { namespace meas { namespace multifit {

class Evaluator : public BaseEvaluator {
public:

    typedef boost::shared_ptr<Evaluator> Ptr;

    Grid::Ptr getGrid() const { return _grid; }
    
    static Ptr make(Grid::Ptr const & grid) {
        return boost::make_shared<Evaluator>(grid);
    }
    static Ptr make(Definition const & definition) {
        return boost::make_shared<Evaluator>(Grid::make(definition));
    }

    virtual double clipToBounds(lsst::ndarray::Array<double,1,1> const & parameters) const {
        return _grid->clipToBounds(parameters.getData());
    }

protected:

    virtual void _evaluateModelMatrix(
        ndarray::Array<double,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & param
    ) const;

    virtual void _evaluateModelMatrixDerivative(
        ndarray::Array<double,3,3> const & modelMatrixDerivative,
        ndarray::Array<double const,2,2> const & modelMatrix,
        ndarray::Array<double const,1,1> const & param
    ) const;

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & param) const;

private:
    
    FRIEND_MAKE_SHARED_1(Evaluator, boost::shared_ptr<lsst::meas::multifit::Grid>);

    explicit Evaluator(Grid::Ptr const & grid);

    Evaluator(Evaluator const & other);
    
    void _initialize();

    Grid::Ptr _grid;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluator
