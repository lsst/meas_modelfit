#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ("fitMixture", "SemiEmpiricalPriorConfig",
           "SoftenedLinearPriorControl")

import numpy as np

from lsst.pex.config import makeConfigClass
from lsst.utils import continueClass

from ..mixture import Mixture
from .priors import (SemiEmpiricalPriorControl, SemiEmpiricalPrior,
                     SoftenedLinearPriorControl, SoftenedLinearPrior,
                     MixturePrior)


SemiEmpiricalPriorConfig = makeConfigClass(SemiEmpiricalPriorControl)

SoftenedLinearPriorConfig = makeConfigClass(SoftenedLinearPriorControl)


@continueClass  # noqa: F811 (FIXME: remove for py 3.8+)
class SemiEmpiricalPrior:  # noqa: F811

    ConfigClass = SemiEmpiricalPriorConfig


@continueClass  # noqa: F811 (FIXME: remove for py 3.8+)
class SoftenedLinearPrior:  # noqa: F811

    ConfigClass = SoftenedLinearPriorConfig


def fitMixture(data, nComponents, minFactor=0.25, maxFactor=4.0,
               nIterations=20, df=float("inf")):
    """Fit a ``Mixture`` distribution to a set of (e1, e2, r) data points,
    returing a ``MixturePrior`` object.

    Parameters
    ----------
    data : numpy.ndarray
        array of data points to fit; shape=(N,3)
    nComponents : int
        number of components in the mixture distribution
    minFactor : float
        ellipticity variance of the smallest component in the initial mixture,
        relative to the measured variance
    maxFactor : float
        ellipticity variance of the largest component in the initial mixture,
        relative to the measured variance
    nIterations : int
        number of expectation-maximization update iterations
    df : float
        number of degrees of freedom for component Student's T distributions
        (inf=Gaussian).
    """
    components = Mixture.ComponentList()
    rMu = data[:, 2].mean()
    rSigma = data[:, 2].var()
    eSigma = 0.5*(data[:, 0].var() + data[:, 1].var())
    mu = np.array([0.0, 0.0, rMu], dtype=float)
    baseSigma = np.array([[eSigma, 0.0, 0.0],
                          [0.0, eSigma, 0.0],
                          [0.0, 0.0, rSigma]])
    for factor in np.linspace(minFactor, maxFactor, nComponents):
        sigma = baseSigma.copy()
        sigma[:2, :2] *= factor
        components.append(Mixture.Component(1.0, mu, sigma))
    mixture = Mixture(3, components, df)
    restriction = MixturePrior.getUpdateRestriction()
    for i in range(nIterations):
        mixture.updateEM(data, restriction)
    return mixture
