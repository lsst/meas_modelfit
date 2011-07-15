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
    
    static Ptr make(Grid::Ptr const & grid, bool usePixelWeights) {
        return boost::make_shared<Evaluator>(grid, usePixelWeights);
    }
    static Ptr make(Definition const & definition, bool usePixelWeights) {
        return boost::make_shared<Evaluator>(Grid::make(definition), usePixelWeights);
    }

    virtual double clipToBounds(lsst::ndarray::Array<double,1,1> const & parameters) const {
        return _grid->clipToBounds(parameters);
    }

    virtual bool checkBounds(lsst::ndarray::Array<double const,1,1> const & parameters) const {
        return _grid->checkBounds(parameters);
    }

protected:

    virtual void _evaluateModelMatrix(
        ndarray::Array<Pixel,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & parameters
    ) const;

#if 0
    virtual void _evaluateModelMatrixDerivative(
        ndarray::Array<Pixel,3,3> const & modelMatrixDerivative,
        ndarray::Array<Pixel const,2,2> const & modelMatrix,
        ndarray::Array<double const,1,1> const & parameters
    ) const;
#endif

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & parameters) const;

private:
    
    FRIEND_MAKE_SHARED_2(Evaluator, boost::shared_ptr<lsst::meas::multifit::Grid>, bool );

    explicit Evaluator(Grid::Ptr const & grid, bool usePixelWeights);
    
    void _initialize();

    bool _usePixelWeights;
    Grid::Ptr _grid;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluator
