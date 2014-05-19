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
#ifndef MEAS_MULTIFIT_Interpreter_h_INCLUDED
#define MEAS_MULTIFIT_Interpreter_h_INCLUDED

#include "lsst/meas/multifit/common.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/Prior.h"

namespace lsst { namespace meas { namespace multifit {

class ModelFitRecord;

class Interpreter {
public:

    PTR(Model) getModel() const { return _model; }

    PTR(Prior) getPrior() const { return _prior; }

    int getParameterDim() const { return _parameterNames.size(); }

    int getNonlinearDim() const { return _model->getNonlinearDim(); }

    int getAmplitudeDim() const { return _model->getAmplitudeDim(); }

    int getFixedDim() const { return _model->getFixedDim(); }

    Model::NameVector const & getParameterNames() const { return _parameterNames; }

    Model::NameVector const & getNonlinearNames() const { return _model->getNonlinearNames(); }

    Model::NameVector const & getAmplitudeNames() const { return _model->getAmplitudeNames(); }

    Model::NameVector const & getFixedNames() const { return _model->getFixedNames(); }

    void packParameters(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & parameters
    ) const;

    void unpackNonlinear(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & nonlinear
    ) const;

    //@{
    /**
     *  @brief Compute quantiles of a single parameter
     *
     *  A quantile is the point at which the cumulative distribution reaches a certain value.
     *  For instance, the median is the quantile at 0.5, and a pair of quantiles at (0.05, 0.95)
     *  is a 90% confidence interval.
     *
     *  @param[in] record          ModelFitRecord from which to extract distributions.
     *  @param[in] fractions       A sorted array of fractions (floating point numbers in the range [0,1])
     *                             at which to compute a quantile value (e.g. 0.5 == median).
     *  @param[in] index           Index of the parameter over which to compute the quantile; an
     *                             integer between 0 and (getNonlinearDim()-1).  All other parameters
     *                             are ignored (meaning they areeffectively marginalized over).
     *
     *  Different functions may be used to operate on the nonlinear parameters, the amplitude parameters,
     *  or the full set of fit parameters.
     */
    virtual ndarray::Array<Scalar,1,1> computeParameterQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions,
        int index
    ) const = 0;

    virtual ndarray::Array<Scalar,1,1> computeNonlinearQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions,
        int index
    ) const = 0;

    virtual ndarray::Array<Scalar,1,1> computeAmplitudeQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions,
        int index
    ) const = 0;
    //@}

    /**
     *  @brief Compute the same quantiles for all nonlinear parameters
     *
     *  Note that this routine computes the quantiles of all N 1-d marginal distributions in turn,
     *  not N-d quantiles (which would need to be contours).
     *
     *  A quantile is the point at which the cumulative distribution reaches a certain value.
     *  For instance, the median is the quantile at 0.5, and a pair of quantiles at (0.05, 0.95)
     *  is a 90% confidence interval.
     *
     *  @param[in] record          ModelFitRecord from which to extract distributions.
     *  @param[in] fractions       A sorted array of fractions (floating point numbers in the range [0,1])
     *                             at which to compute a quantile value (e.g. 0.5 == median).
     *
     *  Different functions may be used to operate on the nonlinear parameters, the amplitude parameters,
     *  or the full set of fit parameters.
     *
     */
    ndarray::Array<Scalar,2,2> computeParameterQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions
    ) const;

    ndarray::Array<Scalar,2,2> computeNonlinearQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions
    ) const;

    ndarray::Array<Scalar,2,2> computeAmplitudeQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions
    ) const;
    //@}

    virtual ndarray::Array<Scalar,1,1> computeParameterMean(ModelFitRecord const & record) const = 0;

    virtual ndarray::Array<Scalar,1,1> computeNonlinearMean(ModelFitRecord const & record) const = 0;

    virtual ndarray::Array<Scalar,1,1> computeAmplitudeMean(ModelFitRecord const & record) const = 0;

    ndarray::Array<Scalar,2,2> computeParameterCovariance(ModelFitRecord const & record) const;

    virtual ndarray::Array<Scalar,2,2> computeParameterCovariance(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & mean
    ) const = 0;

    ndarray::Array<Scalar,2,2> computeNonlinearCovariance(ModelFitRecord const & record) const;

    virtual ndarray::Array<Scalar,2,2> computeNonlinearCovariance(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & mean
    ) const = 0;

    ndarray::Array<Scalar,2,2> computeAmplitudeCovariance(ModelFitRecord const & record) const;

    virtual ndarray::Array<Scalar,2,2> computeAmplitudeCovariance(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & mean
    ) const = 0;

    virtual ~Interpreter() {}

protected:

    Interpreter(Model::NameVector parameterNames, PTR(Model) model, PTR(Prior) prior);

    virtual void _packParameters(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & parameters
    ) const = 0;

    virtual void _unpackNonlinear(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & nonlinear
    ) const = 0;

private:
    Model::NameVector _parameterNames;
    PTR(Model) _model;
    PTR(Prior) _prior;
};

}}} // namespace lsst::afw::table

#endif // !MEAS_MULTIFIT_Interpreter_h_INCLUDED
