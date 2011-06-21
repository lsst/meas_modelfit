from lsst.meas.multifitData import DatasetMapper
import lsst.meas.multifit as mf
from lsst.daf.persistence import ButlerFactory
from lsst.pex.policy import Policy
import lsst.afw.image
import lsst.afw.detection
import lsst.meas.algorithms
import sys

def test():
    bf = ButlerFactory(mapper=DatasetMapper())
    butler = bf.create()

    algorithm = "SHAPELET_MODEL"
    policy = Policy()
    policy.add(algorithm + ".enabled", true)


    results = []
    for i in range(2):
        print >> sys.stderr, "dataset:", i


        if not butler.datasetExists("src", id=i):
            continue
        
        psf = butler.get("psf", id=i)

        exp = butler.get("exp", id=i)
        exp.setPsf(psf)
        sources = butler.get("src", id=i)

        measurePhotometry = lsst.meas.algorithms.makeMeasurePhotometry(exp)
        measurePhotometry.addAlgorithm(algorithm)
        measurePhotometry.configure(policy)
        
        for j, s in enumerate(sources):

            photom = measurePhotometry.measure(lsst.afw.detection.Peak(), s).find(algorithm)
            status = int(photom.get("status"))

            print >> sys.stderr, "source:", j, "flux:", photom.getFlux(), "status:", status


if __name__ == '__main__':
    test()




    
