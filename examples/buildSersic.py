import lsst.meas.multifit
import lsst.afw.geom.ellipses
import lsst.afw.detection
import lsst.afw.math.shapelets
import numpy
from matplotlib import pyplot

# ------- These are the SDSS prescriptions for modified exponential and deVaucouleur profiles ---------

DEFAC = -7.66925
DEVOUT = 8.0
DEVCUT = 7.0
EXPFAC = -1.67835
EXPOUT = 4.0
EXPCUT = 3.0

def deV_profile(r):
    """Truncated de Vaucouleur - copied from SDSS"""
    p = numpy.exp(DEFAC * ((r**2 + 0.0004)**0.125 - 1.0))
    big = r > DEVCUT
    scr = (r[big] - DEVCUT) / (DEVOUT - DEVCUT)
    scr = 1.0 - scr**2
    p[big] *= scr*scr
    p[r > DEVOUT] = 0.0
    return p

def exp_profile(r):
    """Truncated exponential - copied from SDSS"""
    p = numpy.exp(EXPFAC * (r - 1.0))
    big = r > EXPCUT
    scr = (r[big] - EXPCUT) / (EXPOUT - EXPCUT);
    scr = 1.0 - scr**2
    p[big] *= scr * scr
    p[r > EXPOUT] = 0.0
    return p

# ----------------------------------------------------------------------------------------------------

def plotBasis(basis, func):
    pyplot.figure()
    radii = numpy.linspace(0, 4, 500)
    profile = numpy.zeros((radii.size, basis.getSize()), dtype=float)
    basis.evaluateRadialProfile(profile, radii)
    pyplot.subplot(3, 1, 1)
    pyplot.plot(radii, profile[:,0], 'k-')
    pyplot.plot(radii, func(radii), 'r--')
    pyplot.subplot(3, 1, 2)
    pyplot.semilogy(radii, profile[:,0], 'k-')
    pyplot.semilogy(radii, func(radii), 'r--')
    pyplot.xlabel("r / r_e")
    pyplot.ylim(1E-4, func(radii).max())
    pyplot.subplot(3, 1, 3)
    pyplot.bar(range(basis.getMapping().shape[0]), basis.getMapping()[:,0])
    pyplot.axis("off")

def makeExponential():
    components = lsst.meas.multifit.CompoundShapeletBuilder.ComponentVector()
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 1.5))
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 4.0))
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 10.00))
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 20.00))
    sersicRadius = 20.0
    maxRadius = 100.0
    matchRadii = numpy.array([0.0])
    builder = lsst.meas.multifit.CompoundShapeletBuilder.approximate(
        lsst.meas.multifit.ProfileFunction.makeTruncatedExponential(),
        components,
        sersicRadius,
        maxRadius,
        matchRadii
        )
    return builder.build()

def makeDeVaucouleur():
    components = lsst.meas.multifit.CompoundShapeletBuilder.ComponentVector()
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 1.00))
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 3.00))
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 10.00))
    components.push_back(lsst.meas.multifit.ShapeletModelBasis.make(0, 30.00))
    matchRadii = numpy.array([0.0])
    sersicRadius = 25.0
    maxRadius = 100.0
    builder = lsst.meas.multifit.CompoundShapeletBuilder.approximate(
        lsst.meas.multifit.ProfileFunction.makeTruncatedDeVaucouleur(),
        components,
        sersicRadius,
        maxRadius,
        matchRadii
        )
    return builder.build()

def main():
    numpy.set_printoptions(suppress=True, linewidth=180)
    expBasis = makeExponential()
    devBasis = makeDeVaucouleur()
    plotBasis(expBasis, exp_profile)
    plotBasis(devBasis, deV_profile)
    pyplot.show()
    return expBasis, devBasis

if __name__ == "__main__":
    expBasis, devBasis = main()
