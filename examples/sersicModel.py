import lsst.meas.multifit as mf
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.geom.ellipses
import lsst.meas.algorithms
import numpy
import lsst.afw.display.ds9 as ds9
from makeImageStack import makeImageStack

def main():
    flux = 5.0
    centroid = lsst.afw.geom.makePointD(0., 0.)
    axes = lsst.afw.geom.ellipses.Axes(3,1,0)
    logShear = lsst.afw.geom.ellipses.LogShear(axes)
    sersicIndex = 1.0
    model = mf.createSersicModel(flux, centroid, axes, sersicIndex)
    exposureList = makeImageStack(model, 1, centroid[0], centroid[1])

    exposure= exposureList.front()
    ds9.mtv(exposure.getMaskedImage(), frame=0, wcs=exposure.getWcs())

    centroidMeasurer = lsst.meas.algorithms.createMeasureCentroid("SDSS")
    centroid = centroidMeasurer.apply(
        exposure.getMaskedImage().getImage(),
        int(exposure.getWidth()/2 + 1 + exposure.getMaskedImage().getX0()),
        int(exposure.getHeight()/2 + 1 + exposure.getMaskedImage().getY0())
    )

    shapeMeasurer = lsst.meas.algorithms.createMeasureShape("SDSS")
    shape = shapeMeasurer.apply(
        exposure.getMaskedImage(),
        centroid.getX(), 
        centroid.getY()
    )

    moments = lsst.afw.geom.ellipses.Quadrupole(
        shape.getMxx(), shape.getMyy(), shape.getMxy()
    )
    axes = lsst.afw.geom.ellipses.Axes(moments)

    print "Centroid = (%.3f, %.3f)"%(centroid.getX(), centroid.getY())
    print "2nd Moments (xx, xy, yy) = (%.4f, %.4f, %.4f)"%\
        (shape.getMxx(), shape.getMxy(), shape.getMyy())
    print "Ellipse Axes (a, b, theta) = (%.4f, %.4f, %.4f)"%\
        (axes[0], axes[1], axes[2])
    print "0th momeht = %.4f"%shape.getM0()

if __name__== "__main__":
    main()
