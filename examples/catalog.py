import lsst.meas.multifitData 
import lsst.afw.detection
import lsst.afw.image
import lsst.meas.multifit
import lsst.meas.algorithms
import lsst.daf.persistence
import numpy
import sys

FLUX_SCHEMA = lsst.afw.detection.Schema("Flux", 0, lsst.afw.detection.Schema.DOUBLE)
FLUX_ERR_SCHEMA = lsst.afw.detection.Schema("FluxErr", 1, lsst.afw.detection.Schema.DOUBLE)
STATUS_SCHEMA = lsst.afw.detection.Schema("Status", 2, lsst.afw.detection.Schema.INT)
E1_SCHEMA = lsst.afw.detection.Schema("E1", 3, lsst.afw.detection.Schema.DOUBLE)
E2_SCHEMA = lsst.afw.detection.Schema("E2", 4, lsst.afw.detection.Schema.DOUBLE)
RADIUS_SCHEMA = lsst.afw.detection.Schema("R", 5, lsst.afw.detection.Schema.DOUBLE)
COEFF_SCHEMA = lsst.afw.detection.Schema("COEFFICIENTS", 6, lsst.afw.detection.Schema.DOUBLE)

ID = 0
PSF_FLUX = 1
PSF_FLUX_ERR = 2
X = 3
Y = 4
IXX = 5
IYY = 6 
IXY = 7
FIT_STATUS = 8
FIT_FLUX = 9
FIT_FLUX_ERR = 10
FIT_E1 = 11
FIT_E2 = 12
FIT_R = 13
FIT_COEFF_START = 14

def fit(datasets=(0,1,2,3,4,5,6,7,8,9), basisSize=8, fitDelta=True):
    bf = lsst.daf.persistence.ButlerFactory( \
            mapper=lsst.meas.multifitData.DatasetMapper())
    butler = bf.create()
    algorithm = "SHAPELET_MODEL_%d"%basisSize
    policy = lsst.pex.policy.Policy()
    policy.add(algorithm + ".enabled", True)
    policy.add(algorithm + ".fitDeltaFunction", fitDelta)
   
    nVals = 14 + basisSize
    if(fitDelta):
        nVals += 1
    full = None
    for d in datasets:
        if not butler.datasetExists("src", id=d):
            continue
        psf = butler.get("psf", id=d)
        expD = butler.get("exp", id=d)
        miD = expD.getMaskedImage()
        miF = lsst.afw.image.MaskedImageF(
                miD.getImage().convertF(), 
                miD.getMask(), 
                miD.getVariance())

        exposure = lsst.afw.image.ExposureF(miF, expD.getWcs())
        exposure.setPsf(psf)
        sources = butler.get("src", id=d)

        measurePhotometry = lsst.meas.algorithms.makeMeasurePhotometry(exposure)
        measurePhotometry.addAlgorithm(algorithm)
        measurePhotometry.configure(policy)



        partial = numpy.zeros((len(sources), nVals), dtype=float)
        for j, src in enumerate(sources):
            partial[j, ID] = src.getId()
            partial[j, PSF_FLUX] = src.getPsfFlux()
            partial[j, PSF_FLUX_ERR] = src.getPsfFluxErr()
            partial[j, X] = src.getXAstrom()
            partial[j, Y] = src.getYAstrom()
            partial[j, IXX] = src.getIxx()
            partial[j, IYY] = src.getIyy()
            partial[j, IXY] = src.getIxy()

            photom = measurePhotometry.measure(lsst.afw.detection.Peak(), src).find(algorithm)
            partial[j,FIT_STATUS] = photom.get(STATUS_SCHEMA)
            partial[j,FIT_FLUX] = photom.get(FLUX_SCHEMA)
            partial[j,FIT_FLUX_ERR] = photom.get(FLUX_ERR_SCHEMA)
            partial[j,FIT_E1] = photom.get(E1_SCHEMA)
            partial[j,FIT_E2] = photom.get(E2_SCHEMA)
            partial[j,FIT_R] = photom.get(RADIUS_SCHEMA)
            
            for i in range(basisSize):
                partial[j, FIT_COEFF_START + i] = photom.get(i, COEFF_SCHEMA)
        
            if(fitDelta):
                partial[j, FIT_COEFF_START + basisSize] = photom.get(basisSize, COEFF_SCHEMA)

        if(full != None):
            numpy.append(full, partial, 0)
        else:
            full = partial

    return full

def run():
    full = fit();
    if(len(sys.argv) > 1):
        filename = sys.argv[1]
    else filename = "catalog.pickle"
    outfile = open(filename, 'wb')
    pickle.dump(full, outfile)
    outfile.close()

if __name__ == '__main__':
    run()

