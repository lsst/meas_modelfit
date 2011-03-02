# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import math

import lsst.pex.harness.stage as harnessStage

from lsst.pex.logging import Log

import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.pex.exceptions as pexExcept
import lsst.meas.algorithms as measAlg

class SourceToDiaSourceStageParallel(harnessStage.ParallelProcessing):
    """
    Description:
       Glue stage for transforming clipboard objects from SourceSet 
       to DiaSourceSet

    Policy Dictionaty:
    lsst/meas/pipeline/SourceDetectionStageDictionary.paf

    Clipboard Input:
    - SourceSet with key specified by policy attribute inputKey
    - PersistableSourceVector with key "persistable_"+inputKey 
    - CCD-based WCS with key specified by policy attribute ccdWcsKey

    ClipboardOutput:
    - DiaSourceSet with key outputKey.
    - PersistableDiaSourceVector with key "persistable_"+outputKey
    """
    def setup(self):
        self.log = Log(self.log, "SourceToDiaSourceStage - parallel")

        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                                                 "SourceToDiaSourceStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        """
        Converting to DiaSource in the worker process
        """
        self.log.log(Log.INFO, "Executing in process")
       
        self.ccdWcs = clipboard.get(self.policy.get("ccdWcsKey"))
        self.ampBBox = clipboard.get(self.policy.get("ampBBoxKey"))

        dataPolicyList = self.policy.getPolicyArray("data")
        for dataPolicy in dataPolicyList:
            inputKey = dataPolicy.getString("inputKey")
            outputKey = dataPolicy.getString("outputKey")

            sourceSet = clipboard.get(inputKey)
            if sourceSet is None:
                self.log.log(Log.FATAL, "No SourceSet with key " + inputKey)
                continue

            diaSourceSet = afwDet.DiaSourceSet()
            for source in sourceSet:
                diaSource = afwDet.makeDiaSourceFromSource(source)

                (ra, dec, raErr, decErr) = self.raDecWithErrs(
                        diaSource.getXFlux(), diaSource.getYFlux(),
                        diaSource.getXFluxErr(), diaSource.getYFluxErr())
                diaSource.setRaFlux(ra); diaSource.setDecFlux(dec)
                diaSource.setRaFluxErr(raErr); diaSource.setDecFluxErr(decErr)

                (ra, dec, raErr, decErr) = self.raDecWithErrs(
                        diaSource.getXAstrom(), diaSource.getYAstrom(),
                        diaSource.getXAstromErr(), diaSource.getYAstromErr())
                diaSource.setRaAstrom(ra); diaSource.setDecAstrom(dec)
                diaSource.setRaAstromErr(raErr); diaSource.setDecAstromErr(decErr)

                # No errors for XPeak, YPeak
                raDec = self.ccdWcs.pixelToSky(
                    diaSource.getXPeak(), diaSource.getYPeak())
                diaSource.setRaPeak(raDec.getLongitude(afwCoord.DEGREES))
                diaSource.setDecPeak(raDec.getLatitude(afwCoord.DEGREES))

                # Simple RA/decl == Astrom versions
                diaSource.setRa(diaSource.getRaAstrom())
                diaSource.setRaErrForDetection(diaSource.getRaAstromErr())
                diaSource.setDec(diaSource.getDecAstrom())
                diaSource.setDecErrForDetection(diaSource.getDecAstromErr())

                diaSourceSet.append(diaSource)

            persistableSet = afwDet.PersistableDiaSourceVector(diaSourceSet)

            clipboard.put(outputKey, diaSourceSet)
            clipboard.put("persistable_" + outputKey, persistableSet)

    def raDecWithErrs(self, x, y, xErr, yErr, pixToSkyAffineTransform=None):
        """Use wcs to transform pixel coordinates x, y and their errors to 
        sky coordinates ra, dec with errors. If the caller does not provide an
        affine approximation to the pixel->sky WCS transform, an approximation
        is automatically computed (and used to propagate errors). For sources
        from exposures far from the poles, a single approximation can be reused
        without introducing much error.

        Note that the affine transform is expected to take inputs in units of
        pixels to outputs in units of degrees. This is an artifact of WCSLIB
        using degrees as its internal angular unit.

        Sky coordinates and their errors are returned in units of degrees.
        """
        ampX = x - self.ampBBox.getX0()
        ampY = y - self.ampBBox.getY0()
        sky = self.ccdWcs.pixelToSky(ampX, ampY)
        if pixToSkyAffineTransform is None:
            pixToSkyAffineTransform = self.ccdWcs.linearizePixelToSky(sky)
        raErr, decErr = self.raDecErrs(xErr, yErr, pixToSkyAffineTransform)
        return (sky.getLongitude(afwCoord.DEGREES),
                sky.getLatitude(afwCoord.DEGREES),
                raErr,
                decErr)

    def raDecErrs(self, xErr, yErr, pixToSkyAffineTransform):
        """Propagates errors in pixel space to errors in ra, dec space
        using an affine approximation to the pixel->sky WCS transform
        (e.g. as returned by lsst.afw.image.Wcs.linearizeAt).

        Note that pixToSkyAffineTransform is expected to take inputs in units
        of pixels to outputs in units of degrees. This is an artifact of WCSLIB
        using degrees as its internal angular unit.

        Errors are returned in units of degrees.
        """
        t = pixToSkyAffineTransform
        varRa  = t[0]**2 * xErr**2 + t[2]**2 * yErr**2
        varDec = t[1]**2 * xErr**2 + t[3]**2 * yErr**2
        return (math.sqrt(varRa), math.sqrt(varDec))


class SourceToDiaSourceStage(harnessStage.Stage):
    parallelClass = SourceToDiaSourceStageParallel
