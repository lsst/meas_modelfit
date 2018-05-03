# This file is part of meas modelfit.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np


def makeGaussian(x, y, scale, muX, muY, varX, varXY, varY):
    """Create an elliptical Gaussian.

    Parameters
    ----------
    x : 2D numpy array
        An array containing the x coordinates the Gaussian will be evaluated on.
        most likely the result of a numpy.indices call
    y : 2D numpy array
        An array containing the y coordinates the Gaussian will be evaluated on.
        most likely the result of a numpy.indices call
    scale : `float`
        The value the resulting Gaussian will have when summed over all pixels.
    muX : `float`
        The central position of the Gaussian in the x direction
    muY : `float`
        The central position of the Gaussian in the y direction
    varX : `float`
        The variance of the Gaussian about the muX position
    varXY : `float`
        The covariance of the Gaussian in x and y
    varY : `float`
        The variance of the Gaussian about the muY position

    Returns
    -------
    Gaussian : 2D numpy array
        The Gaussian array generated from the input values
    """

    rho = varXY/(varX**0.5*varY**0.5)

    psf = np.exp(-1/(2*(1-rho**2)) *
                 ((x-muX)**2/varX+(y - muY)**2/varY -
                 2*rho*(x-muX)*(y-muY)/(varX**0.5*varY**0.5)))

    psf /= psf.sum()
    return scale*psf


def buildUncertanty(imShape, W, uncertanty):
    """ Propagate pixel uncertainties to uncertainties in weighted moments

    Parameters
    ----------
    imShape : tuple(float, float)
        The shape of image for which weighted moments have been calculated
    W : iterable
        An iterable object with six elements corresponding to the moments used
        in the weighted moment calculation, scale, mean in x, mean in y, variance
        in x, covariance between x and y, and variance in y.
    uncertanty : `float`
        Uncertainty in the pixel value. This is a single value, as this routine
        assumes errors are background dominated, and uncorrelated

    Returns
    -------
    covarianceMatrix : 2D 6x6 numpy array of floats
        This is the covariance matrix on the measured moments with uncertainties
        propagated from pixel uncertainties
    """
    yInd, xInd = np.indices(imShape)
    weightImage = makeGaussian(xInd, yInd, *W)
    sigmaImage = np.eye(weightImage.size)*uncertanty
    MomentWeightMatrix = np.zeros((6, weightImage.size))

    weightImageFlat = weightImage.flatten()
    xIndFlat = xInd.flatten()
    yIndFlat = yInd.flatten()

    MomentWeightMatrix[0] = weightImageFlat
    MomentWeightMatrix[1] = weightImageFlat*xIndFlat
    MomentWeightMatrix[2] = weightImageFlat*yIndFlat
    MomentWeightMatrix[3] = weightImageFlat*xIndFlat**2
    MomentWeightMatrix[4] = weightImageFlat*xIndFlat*yIndFlat
    MomentWeightMatrix[5] = weightImageFlat*yIndFlat**2

    return np.dot(MomentWeightMatrix, np.dot(sigmaImage, np.transpose(MomentWeightMatrix)))


def measureMoments(image, W):
    """ Calculate weighted moments of the input image with the given weight array

    Parameters
    ----------
    image : 2D numpy array of floats
        This is the input postage stamp of a source for which the weighted moments are
        to be measured
    W : 2D numpy array of floats
        Array of floats that are used as weights when calculating moments on the input
        image. Array must be the same shape image

    Returns
    -------
    moments : 6 element numpy array
        These are the weighted moments as measured from the input image in the order of
        0th, 1st X, 1st Y, 2nd X, 2nd XY, 2nd Y

    Raises
    ------
    AssertionError: Raises if the input arrays are not the same shape
    """
    assert image.shape == W.shape, "Input image and weight array must be the same shape"

    yInd, xInd = np.indices(image.shape)
    weightImage = makeGaussian(xInd, yInd, *W)

    zero = np.sum(image*weightImage)
    oneX = np.sum(image*weightImage*xInd)
    oneY = np.sum(image*weightImage*yInd)
    twoX = np.sum(image*weightImage*xInd**2)
    twoXY = np.sum(image*weightImage*xInd*yInd)
    twoY = np.sum(image*weightImage*yInd**2)

    return np.array((zero, oneX, oneY, twoX, twoXY, twoY)), weightImage
