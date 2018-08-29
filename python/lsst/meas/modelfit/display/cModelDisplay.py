__all__ = ("displayReconstructedCmodel", "displayReconstrucedCmodelMpl",
           "buildCModelImages", "reconstructCModel")

import numpy as np
import matplotlib.pyplot as plt

import lsst.afw.image as afwImage
import lsst.meas.modelfit as measMod
import lsst.shapelet as shapelet


def displayReconstructedCmodel(exposure, record, config, display="mpl"):
    """ Display an image, the Cmodel, and residuals

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object that contains the source which was modeled
    record : `lsst.afw.table.SourceRecord`
        Record object which contains the measurements made on the source
    config : `lsst.meas.modelfit.CModel(SingleFrame/Forced)Config`
        Configuration object of the CModel plugin used in the measurement
        process
    display : `str`, optional
        Display in which to render the data, may be mpl for matplotlib or afw
        for afwDisplay, defaults to mpl

    Raises
    ------
    `NotImplementedError`
        If the display backend specified is not implemented
    """

    if display is "mpl":
        displayReconstrucedCmodelMpl(exposure, record, config)
    elif display is "afw":
        raise NotImplementedError("The afw display backend is not yet "
                                  "implemented")
    else:
        raise NotImplementedError("The display backend '{}' is not a supported"
                                  "backend".format(display))


def displayReconstrucedCmodelMpl(exposure, record, config):
    """ Display an image, the Cmodel, and residuals using matplotlib

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object that contains the source which was  modeled
    record : `lsst.afw.table.SourceRecord`
        Record object which contains the measurements made on the source
    config : `lsst.meas.modelfit.CModel(SingleFrame/Forced)Config`
        Configuration object of the CModel plugin used in the measurement
        process
    """

    subImage, devIm, expIm, jointIm = buildCModelImages(exposure, record,
                                                        config)
    # Get the min and max values for the sub image to use for scaling
    subImMin = subImage.array.min()
    subImMax = subImage.array.max()

    # Calculate the scaling to use for the residual images
    devResidual = subImage.array-devIm.array
    expResidual = subImage.array-expIm.array
    jointResidual = subImage.array-jointIm.array
    residualList = (devResidual, expResidual, jointResidual)

    differences = [(x.max()-x.max()) for x in residualList]
    maxDifferenceIndex = np.argmax(differences)
    resImMin = residualList[maxDifferenceIndex].min()
    resImMax = residualList[maxDifferenceIndex].max()

    # Build the image figure
    fig, axes = plt.subplots(3, 3, sharex='col', sharey='row', figsize=(8, 5))
    fig.subplots_adjust(left=0.25, right=0.8)
    lCBar = fig.add_axes([0.1, 0.15, 0.05, 0.7])
    rCBar = fig.add_axes([0.85, 0.15, 0.05, 0.7])

    # Populate just the exposures in the appropriate places
    for i in range(3):
        axes[i, 0].imshow(subImage.array, vmin=subImMin, vmax=subImMax)

    # Populate dev panels
    axes[0, 1].imshow(devIm.array, vmin=subImMin, vmax=subImMax)
    axes[0, 2].imshow(devResidual, vmin=resImMin, vmax=resImMax, cmap="BrBG")

    # Populate exp panels
    axes[1, 1].imshow(expIm.array, vmin=subImMin, vmax=subImMax)
    axes[1, 2].imshow(expResidual, vmin=resImMin, vmax=resImMax, cmap="BrBG")

    # Populate joint panels
    axes[2, 1].imshow(jointIm.array, vmin=subImMin, vmax=subImMax)
    axes[2, 2].imshow(jointResidual, vmin=resImMin, vmax=resImMax, cmap="BrBG")

    axes[0, 0].set_title("Image")
    axes[0, 1].set_title("Model")
    axes[0, 2].set_title('Residuals')

    axes[0, 0].set_ylabel("Dev")
    axes[1, 0].set_ylabel("Exp")
    axes[2, 0].set_ylabel("Joint")

    fig.colorbar(axes[0, 0].get_images()[0], lCBar)
    lCBar.yaxis.set_ticks_position('left')
    fig.colorbar(axes[maxDifferenceIndex, 2].get_images()[0], rCBar)

    plt.show()


def buildCModelImages(exposure, record, config):
    """ Create Images out of the CModel for the given record

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object that contains the source which was  modeled
    record : `lsst.afw.table.SourceRecord`
        Record object which contains the measurements made on the source
    config : `lsst.meas.modelfit.CModel(SingleFrame/Forced)Config`
        Configuration object of the CModel plugin used in the measurement
        process

    Returns
    -------
    subImage : `lsst.afw.image.ImageF`
        Sub image of the original data taken from a region defined by the
        bounding box of the footprint for the object defined in the given
        source record
    devIm : `lsst.afw.image.ImageF`
        Image created from the dev component of the CModel for the supplied
        record at the same pixel locations as subImage
    expIm: `lsst.afw.image.ImageF`
        Image created from the exp component of the CModel for the supplied
        record at the same pixel locations as subImage
    jointIm :
        Image created from the joint fit of the dev and exp components of the
        CModel for the supplied record at the same pixel locations as subImage
    """

    dev, exp, jointDev, jointExp = reconstructCModel(exposure, record, config)
    # Get exposure cutout
    footBBox = record.getFootprint().getBBox()
    subImage = afwImage.ImageF(exposure.getImage(), footBBox)

    # Build the psf
    shapeletPsfKey = shapelet.MultiShapeletFunctionKey(
        record.schema[config.psfName])
    psfApprox = record.get(shapeletPsfKey)

    # Build the dev Image from the shapelet function
    devIm = afwImage.ImageF(footBBox)
    dev = dev.convolve(psfApprox)
    dev.evaluate().addToImage(devIm)

    # Build the exp image from the shapelet function
    expIm = afwImage.ImageF(footBBox)
    exp = exp.convolve(psfApprox)
    exp.evaluate().addToImage(expIm)

    # Build the joint image from the shapelet function
    jointIm = afwImage.ImageF(footBBox)
    jointDev = jointDev.convolve(psfApprox)
    jointExp = jointExp.convolve(psfApprox)
    jointDev.evaluate().addToImage(jointIm)
    jointExp.evaluate().addToImage(jointIm)

    return subImage, devIm, expIm, jointIm


def reconstructCModel(exposure, record, config):
    """ Reconstruct the CModel for the given record

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object that contains the source which was  modeled
    record : `lsst.afw.table.SourceRecord`
        Record object which contains the measurements made on the source
    config : `lsst.meas.modelfit.CModel(SingleFrame/Forced)Config`
        Configuration object of the CModel plugin used in the measurement
        process

    Returns
    -------
    devShapelet : `lsst.shapelet.MultiShapeletFunction`
        Multi-component shapelet model of the dev component of CModel
    expShapelet : `lsst.shapelet.MultiShapeletFunction`
        Multi-component shapelet model fo the exp component of CModel
    devJointShapelet : `lsst.shapelet.MultiShapeletFunction`
        Multi-component shapelet model of the dev component of CModel jointly
        fit with the exp component
    expJointShapelet : `lsst.shapelet.MultiShapeletFunction`
        Multi-component shapelet model of the exp component of Cmodel jointly
        fit with the dev component
    """

    # build a unit system transformation object
    center = record.getCentroid()
    position = exposure.getWcs().pixelToSky(center)
    measSys = measMod.UnitSystem(exposure)
    approxFlux = record.get("base_PsfFlux_instFlux")
    fitSys = measMod.UnitSystem(position, exposure.getCalib(), approxFlux)
    fitSysToMeasSys = measMod.LocalUnitTransform(center, fitSys, measSys)

    # Build the Shapelet objects
    ctrl = config.makeControl()
    baseName = "modelfit_CModel"
    nonlinearKeys = ["{}_{{model}}_nonlinear_{p}".format(baseName, p=p)
                     for p in range(3)]
    fixedKeys = ["{}_{{model}}_fixed_{p}".format(baseName, p=p)
                 for p in range(2)]
    fluxKey = "{}_{{model}}_instFlux".format(baseName)

    # fetch the aperture corrections, if this fails set it to one
    try:
        apCorr = record.get("{}_apCorr".format(baseName))
    except Exception:
        print("Warning, problem retrieving aperture correction, using a value"
              " of 1")
        apCorr = 1

    # Get parameters for the dev model
    devNonLinearParams = np.array([record.get(key.format(model="dev"))
                                  for key in nonlinearKeys])
    devFixedParams = np.array([record.get(key.format(model="dev"))
                              for key in fixedKeys])
    devAmp = np.array([record.get(fluxKey.format(model="dev"))])
    devAmp /= apCorr
    devShapelet = ctrl.dev.getModel().makeShapeletFunction(devNonLinearParams,
                                                           devAmp,
                                                           devFixedParams)
    devShapelet.transformInPlace(fitSysToMeasSys.geometric)

    # Get parameters for the exp model
    expNonLinearParams = np.array([record.get(key.format(model="exp"))
                                  for key in nonlinearKeys])
    expFixedParams = np.array([record.get(key.format(model="exp"))
                              for key in fixedKeys])
    expAmp = np.array([record.get(fluxKey.format(model="exp"))])
    expAmp /= apCorr
    expShapelet = ctrl.exp.getModel().makeShapeletFunction(expNonLinearParams,
                                                           expAmp,
                                                           expFixedParams)
    expShapelet.transformInPlace(fitSysToMeasSys.geometric)

    # Get joint shapelet model
    fracDev = record.get("{}_fracDev".format(baseName))
    jointFlux = np.array([record.get("{}_instFlux".format(baseName))])
    jointFlux /= apCorr
    devJointShapelet = ctrl.dev.getModel()\
        .makeShapeletFunction(devNonLinearParams, jointFlux*fracDev,
                              devFixedParams)
    devJointShapelet.transformInPlace(fitSysToMeasSys.geometric)

    expJointShapelet = ctrl.exp.getModel()\
        .makeShapeletFunction(expNonLinearParams, jointFlux*(1-fracDev),
                              expFixedParams)
    expJointShapelet.transformInPlace(fitSysToMeasSys.geometric)

    return devShapelet, expShapelet, devJointShapelet, expJointShapelet
