from lsst.meas.multifitData import DatasetMapper
import lsst.meas.multifit as mf
from lsst.daf.persistence import ButlerFactory
from lsst.pex.policy import Policy

def test():
    bf = ButlerFactory(mapper=DatasetMapper())
    butler = bf.create()

    policy = Policy()
    policy.add("nGrowFp", 3)
    policy.add("isVariable", False)
    policy.add("isPositionActive", False)
    policy.add("isRadiusActive", True)
    policy.add("isEllipticityActive", True)
    policy.add("maskPlaneName", "BAD")
    policy.add("basisName", "ed+15:4000")

    results = []
    for i in range(2):
        if not butler.datasetExists("src", id=i):
            continue
        psf = butler.get("psf", id=i)
        exposure = butler.get("exp", id=i)
        exposure.setPsf(psf)
        sources = butler.get("src", id=i)
        results.append(mf.utils.processExposure(exposure, sources, policy))

if __name__ == '__main__':
    test()




    
