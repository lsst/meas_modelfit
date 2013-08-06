#!/usr/bin/env python

import sys

import lsst.afw.image
import lsst.afw.geom
import lsst.daf.persistence
import lsst.meas.multifit

dataId = dict(visit=100, raft="2,2", sensor="1,1")

def main(root, indices=(0,)):
    butler = lsst.daf.persistence.Butler(root)
    dataRef = butler.dataRef("calexp", dataId=dataId)
    inExp = dataRef.get("calexp", immediate=True)
    inCat = dataRef.get("modelfits", immediate=True)
    outCat = lsst.meas.multifit.ModelFitCatalog(inCat.getTable())
    bbox = lsst.afw.geom.Box2I()
    for index in indices:
        outCat.append(inCat[index])
        bbox.include(inCat[index].getFootprint().getBBox())
    bbox.grow(2)
    bbox.clip(inExp.getBBox(lsst.afw.image.PARENT))
    outExp = inExp.Factory(inExp, bbox, lsst.afw.image.PARENT, True)
    outExp.writeFits("calexp.fits")
    outCat.writeFits("outcat.fits")

if __name__ == "__main__":
    main(root=sys.argv[1], indices=tuple(int(s) for s in sys.argv[2:]))
