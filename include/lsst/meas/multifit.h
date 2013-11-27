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

#ifndef LSST_MEAS_MULTIFIT_H
#define LSST_MEAS_MULTIFIT_H

#include "lsst/meas/multifit/MarginalSampling.h"
#include "lsst/meas/multifit/DirectSampling.h"
#include "lsst/meas/multifit/ModelFitRecord.h"
#include "lsst/meas/multifit/AdaptiveImportanceSampler.h"
#include "lsst/meas/multifit/Sampling.h"
#include "lsst/meas/multifit/TruncatedGaussian.h"
#include "lsst/meas/multifit/Likelihood.h"
#include "lsst/meas/multifit/ProjectedLikelihood.h"
#include "lsst/meas/multifit/UnitSystem.h"
#include "lsst/meas/multifit/Interpreter.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/integrals.h"
#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/Mixture.h"
#include "lsst/meas/multifit/optimizer.h"
#include "lsst/meas/multifit/psf.h"

#endif // !LSST_MEAS_MULTIFIT_H
