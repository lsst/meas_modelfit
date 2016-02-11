#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import glob
import unittest
import math
import numpy
import os
import pyfits
import time

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.afw.table
import lsst.afw.detection
import lsst.meas.modelfit
import lsst.meas.base
import lsst.meas.algorithms



def makePsf(data, max=None):
    if not max is None:
        trim0 = (data.shape[0] - max)/2
        trim1 = (data.shape[1] - max)/2
        if trim0 < 0 or trim1 < 0:
            raise BaseException("psfSize is smaller than the psf supplied")
        data = data[trim0:trim0+max, trim1:trim1+max]
    kernel = lsst.afw.math.FixedKernel(lsst.afw.image.ImageD(data))
    return lsst.meas.algorithms.KernelPsf(kernel)


def getStats(self, vals, thresh=None, iters=1):
    clip = None
    clipped = []
    for i in range(iters):
        if i == 0: clip = None
        else: clip = valAvg + thresh * valStdev
        valSum = 0.0
        valSS = 0.0
        valCount= 0
        for val in vals: 
            if clip is None or val < clip:
                valSum = valSum + val 
                valSS = valSS + val * val
                valCount = valCount + 1
            else:
                if i == (iters-1):
                    clipped.append(val)
        valAvg = valSum/valCount
        valStdev = math.sqrt((valCount * valSS - valSum * valSum)/ (valCount * (valCount - 1)))
        if thresh == None:
             return valAvg, valStdev
    return valAvg, valStdev, clip, clipped

def runTimeTest(psftype, libname, model, exposure, psfSize=None, indices=None, verbose=False, timelimit=None, start=0, count=10000):
    config = lsst.meas.base.SingleFrameMeasurementTask.ConfigClass()
    config.slots.centroid = "base_GaussianCentroid"
    config.slots.shape = None
    config.slots.psfFlux = None
    config.slots.apFlux = None
    config.slots.instFlux = None
    config.slots.modelFlux = None
    config.slots.calibFlux = None
    config.doReplaceWithNoise = False
    config.plugins.names = ["base_GaussianCentroid", "base_PsfFlux", "base_SdssShape", "modelfit_ShapeletPsfApprox", "modelfit_CModel"]
    config.plugins["modelfit_CModel"].psfName=model
    config.plugins["modelfit_ShapeletPsfApprox"].sequence = [model]
    config.plugins["modelfit_ShapeletPsfApprox"].models[model].optimizer.maxOuterIterations = 2000
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=schema)
    task.log.setThreshold(task.log.WARN)
    measCat = lsst.afw.table.SourceCatalog(schema)
    measRecord = measCat.addNew()
    flags = ("flag", "flag_max_inner_iterations", "flag_max_outer_iterations", "flag_contains_nan", "flag_exception")
    timesSPA = []
    flagCounts = {} 
    excepts = []
    flagKeys = {}
    processed = 0 
    nGood = 0
    Sum = 0
    SumSq = 0
    SumCM = 0
    SumSqCM = 0
    for i, flag in enumerate(flags):
        flagCounts[flag] = []
        flagKeys[flag] = measRecord.getSchema().find("modelfit_ShapeletPsfApprox_%s_%s"%(model, flag)).key
    hdus = pyfits.open("data/%s/psf_library_%s.fits"%(psftype,libname))

    # create the output directory
    outdir = "results/%s/%s"%(psftype, model)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print "xxmodel\tpsftype\tlib\tpsfSize\tn\tnGood\tSPAavg\tSPAstd\tCMavg\tCMstd"
    for i in range(len(hdus)):
        if i < start:
            continue
        if i >= count:
            break
        if not indices is None and not i in indices:
            continue
        exposure.setXY0(lsst.afw.geom.Point2I(0, 0))
        psfData = hdus[i].data.astype(numpy.float64)
        psf = makePsf(psfData, max=psfSize)
        psfWidth = psfSize
        psfHeight = psfSize
        exposure.setPsf(psf)

        CD = numpy.array([[5.55E-5, 0.0], [0.0, 5.55E-5]])
        crpix = lsst.afw.geom.Point2D(0.0,0.0)
        crval = lsst.afw.geom.Point2D(0.0,0.0)
        exposure.setWcs(lsst.afw.image.Wcs(crval, crpix, CD))

        exposure.setXY0(lsst.afw.geom.Point2I(0, 0))
        dettask = lsst.meas.algorithms.SourceDetectionTask()
        dettask.log.setThreshold(dettask.log.WARN)
        footprints = dettask.detectFootprints(exposure, sigma=4.0).positive.getFootprints()
        measRecord.setFootprint(footprints[0])
        for plugin in task.plugins.keys():
            startTime = time.time()
            try:
                task.plugins[plugin].measure(measRecord, exposure)
            except:
                excepts.append(i)
                pass

            if plugin == "modelfit_ShapeletPsfApprox":
                timeSPA = time.time() - startTime
                lastFlag = ""
                flagSet = False
                for flag in flags:
                    if measRecord.get(flagKeys[flag]):
                        flagSet = True
                        flagCounts[flag].append(i)
                        measRecord.set(flagKeys[flag], False)
                        lastFlag = flag
                if not flagSet:
                    timesSPA.append(timeSPA)
                    nGood = nGood + 1
                    Sum = Sum + timeSPA
                    SumSq = SumSq + (timeSPA * timeSPA)
            if plugin == "modelfit_CModel":
                timeCM = time.time() - startTime
                if not flagSet:
                    SumCM = SumCM + timeCM
                    SumSqCM = SumSqCM + (timeCM * timeCM)

        if verbose:
            print "%d of %d"%(i, len(hdus)), " SPA=", timeSPA, " CM=", timeCM, " %s"%lastFlag
        processed = processed + 1

    Avg = Sum/nGood
    AvgCM = SumCM/nGood
    print "%.3f"%((SumSq/nGood) - (Avg*Avg)), "%.3f"%(AvgCM), "%.3f"%((SumSqCM/nGood) - (AvgCM*AvgCM))
    print "xx", model, psftype, libname, psfSize, nGood, processed, "%.3f"%(Avg), "%.3f"%(math.sqrt((SumSq/nGood) - (Avg*Avg))), "%.3f"%(AvgCM), "%.3f"%(math.sqrt((SumSqCM/nGood) - (AvgCM*AvgCM)))
    for flag in flags:
        if len(flagCounts[flag]) > 0:
            print "Failure counts for ", flag, " = ", flagCounts[flag]
    return 

if __name__ == "__main__":
    for model in ("Full", ):
        exposure = lsst.afw.image.ExposureF(os.path.join("data", "exposure.fits"))
        for libname in ("12",):
          for seeing in ("0.5",):
            for size in (25, 33):
              runTimeTest("f2_%s"%seeing, libname, model, exposure, count=50, psfSize=size, verbose=True)
