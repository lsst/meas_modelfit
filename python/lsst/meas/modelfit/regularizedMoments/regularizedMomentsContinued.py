import numpy as np


def makeGaussian(x, y, scale, muX, muY, varX, varXY, varY):
    rho = varXY/(varX**0.5*varY**0.5)
    norm = 1/(2*np.pi*varX*varY*(1-rho**2)**0.5)

    psf = np.exp(-1/(2*(1-rho**2)) *
                      ((x-muX)**2/varX+(y - muY)**2/varY -
                      2*rho*(x-muX)*(y-muY)/(varX**0.5*varY**0.5)))
    
    psf /= psf.sum()
    psf = np.zeros(x.shape)
    for i in range(y.shape[0]):
        for j in range(x.shape[1]):
            mu = np.array([x[i, j] - muX, y[i, j] - muY])
            sig = np.array([[varX, varXY], [varXY, varY]])
            temp = np.dot(np.linalg.inv(sig), mu)
            psf[i, j] = np.exp(-1/2.*np.dot(mu, temp))
    psf /= psf.sum()
    return scale*psf


def buildUncertanty(imShape, W, uncertanty):
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
    yInd, xInd = np.indices(image.shape)
    weightImage = makeGaussian(xInd, yInd, *W)

    zero = np.sum(image*weightImage)
    oneX = np.sum(image*weightImage*xInd)
    oneY = np.sum(image*weightImage*yInd)
    twoX = np.sum(image*weightImage*xInd**2)
    twoXY = np.sum(image*weightImage*xInd*yInd)
    twoY = np.sum(image*weightImage*yInd**2)

    return np.array((zero, oneX, oneY, twoX, twoXY, twoY)), weightImage
