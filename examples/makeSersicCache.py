import lsst.afw.image as afwImage
import lsst.meas.multifit as measMult
import lsst.pex.policy as pexPol
import lsst.afw.display.ds9 as ds9
import sys

def makeCache():
    pol = pexPol.Policy()
    pol.set("kMax", 10.0)
    pol.set("kResolution", 1.0)
    cache = measMult.makeSersicCache(pol)

    data = cache.getDataPoints()
    image = afwImage.ImageD(data.shape[0], data.shape[1])
    for i in range(image.getWidth()):
        for j in range(image.getHeight()):
            image.set(i, j, data[i,j])

    ds9.mtv(image, frame=0)

if __name__ == '__main__':
    makeCache()
