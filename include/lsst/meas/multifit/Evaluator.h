// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#include "multifit/Grid.hpp"
#include "multifit/Interpreter.hpp"
#include "multifit/ParameterConstraint.hpp"
#include <modeling/parameters.hpp>
#include <modeling/Evaluator.hpp>

namespace lsst { namespace meas { namespace multifit {

class Evaluator : public modeling::Evaluator {
public:

    typedef boost::shared_ptr<Evaluator> Ptr;

    boost::shared_ptr<agl::wcs::Projection const> getWCS() const;

    Definition makeDefinition() const;
    Definition makeDefinition(parameters::ConstArray const & parameters) const;
    
    Interpreter makeInterpreter(std::string const & ellipse_name) const {
        return Interpreter(_grid, ellipse_name);
    }

    ParameterConstraint::Ptr makeParameterConstraint(double max_centroid_shift) const {
        return boost::make_shared<ParameterConstraint>(max_centroid_shift, _grid);
    }

    static Ptr make(Definition const & definition);

protected:

    virtual void _evaluateModelMatrix(
        ndarray::Array<double,2,2> const & matrix,
        modeling::parameters::ConstArray const & param
    ) const;

    virtual void _evaluateModelDerivative(
        ndarray::Array<double,3,3> const & matrix,
        modeling::parameters::ConstArray const & param
    ) const;

    virtual void _writeInitialParameters(parameters::Array const & param) const;

private:
    
    FRIEND_MAKE_SHARED_1(multifit::Evaluator, boost::shared_ptr<multifit::Grid>);

    explicit Evaluator(boost::shared_ptr<Grid> const & grid);

    Evaluator(Evaluator const & other);
    
    void _initialize();

    boost::shared_ptr<Grid> _grid;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluator
