from lsst.meas.multifitData import DatasetMapper
import lsst.meas.multifit as mf
from lsst.daf.persistence import ButlerFactory
from lsst.pex.policy import Policy

def test():
    bf = ButlerFactory(mapper=DatasetMapper())
    butler = bf.create()

    policy = Policy()
    policy.add("isVariable", False)
    policy.add("isPositionActive", False)
    policy.add("isRadiusActive", False)
    policy.add("isEllipticityActive", False)
    policy.add("maskPlaneName", "BAD")


    results = {}
    dsTypes = ["highPs", "lowPs", "highSg", "lowSg"]    
    for i in range(20):
        psf = butler.get("psf", id=i)
        for ds in dsTypes:
            if not butler.datasetExists("src", id=i, dsType=ds):
                continue

            src = butler.get("src", id=i, dsType=ds)
            cutout = butler.get("exp", id=i, dsType=ds)
            cutout.setPsf(psf)

            result = mf.utils.fitSource(cutout, src, policy)
            results[ds+i]=result


if __name__ == '__main__':
    test()




    
