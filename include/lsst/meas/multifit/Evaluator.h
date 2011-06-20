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

    /// @brief Size of data vector (number of rows of model matrix).
    virtual int getPixelCount() const { return _grid->getPixelCount(); }

    /// @brief Size of coefficient vector (number of colums of model matrix).
    virtual int getCoefficientCount() const { return _grid->getCoefficientCount(); }

    /// @brief Number of parameters.
    virtual int getParameterCount() const { return _grid->getParameterCount(); }

    /**
     *  @brief Return the natural log of the normalization term of the Gaussian likelihood.
     *
     *  This is a constant that does not depend on the parameters or coefficients, and only
     *  matters when the Bayesian evidence is computed.  For per-pixel uncertainties @f$\sigma_i@f$
     *  or, equivalently, a diagonal covariance matrix @f$\Sigma@f$, the returned value is
     *  @f[
     *    \sum_i \ln \frac{2\pi}{\sigma_i} = \frac{1}{2}\ln \left|2\pi\Sigma^{-1}\right|
     *  @f]
     */
    virtual double getLogPixelErrorSum() const { return _logPixelErrorSum; }

    /**
     *  @brief Data vector.
     *
     *  If the data vector is weighted (divided by sigma) the evaluted model matrix should be as well.
     */
    virtual lsst::ndarray::Array<Pixel const,1,1> getDataVector() const { return _dataVector; }

    Grid::Ptr getGrid() const { return _grid; }
    
    static Ptr make(Grid::Ptr const & grid) {
        return Ptr(new Evaluator(grid));
    }
    static Ptr make(Definition const & definition) {
        return Ptr(new Evaluator(Grid::make(definition)));
    }

protected:

    virtual double _clipToBounds(ndarray::Array<double,1,1> const & parameters) const {
        return _grid->clipToBounds(parameters);
    }

    virtual bool _checkBounds(lsst::ndarray::Array<double const,1,1> const & parameters) const {
        return _grid->checkBounds(parameters);
    }

    virtual CoefficientPrior::ConstPtr _evaluate(
        ndarray::Array<Pixel,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & parameters
    ) const;

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & parameters) const;

private:

    explicit Evaluator(Grid::Ptr const & grid);

    Evaluator(Evaluator const & other);

    Grid::Ptr _grid;
    double _logPixelErrorSum;
    ndarray::Array<Pixel,1,1> _dataVector;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluator
