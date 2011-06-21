#!/usr/bin/env python

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
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as measAlgorithms
import lsst.meas.multifit as mf
import lsst.pex.policy as policy
import lsst.utils.tests as utilTests
from lsst.daf.persistence import ButlerFactory
import unittest
import eups
import sys
import math
import numpy

if(eups.productDir("meas_multifitData")):
    import lsst.meas.multifitData


class GaussianPsfTestCase(unittest.TestCase):
    """A test case detecting and measuring Gaussian PSFs"""
    def setUp(self):
        FWHM = 5
        psf = afwDetection.createPsf("SingleGaussian", 15, 15, FWHM/(2*math.sqrt(2*math.log(2))))
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 100))

        self.xc, self.yc, self.flux = 45, 55, 1000.0
        mi.getImage().set(self.xc, self.yc, self.flux)
        mi.getVariance().set(1e-6)

        cnvImage = mi.Factory(mi.getDimensions())
        afwMath.convolve(cnvImage, mi, psf.getKernel(), afwMath.ConvolutionControl())
        self.exp = afwImage.makeExposure(cnvImage)
        self.exp.setPsf(psf)

        axes = afwGeom.ellipses.Axes(1, 1, 0)
        quad = afwGeom.ellipses.Quadrupole(axes)
        point = afwGeom.Point2D(self.xc, self.yc)
        ellipse = afwGeom.ellipses.Ellipse(axes, point)        
        print >>sys.stderr, ellipse
        self.source = afwDetection.Source()
        self.source.setXAstrom(self.xc)
        self.source.setYAstrom(self.yc)
        self.source.setIxx(quad.getIXX())
        self.source.setIyy(quad.getIYY())
        self.source.setIxy(quad.getIXY())
        self.source.setFootprint(afwDetection.Footprint(ellipse))

        self.mp = measAlgorithms.makeMeasurePhotometry(self.exp)
        self.mp.addAlgorithm("SHAPELET_MODEL")

        self.pol = policy.Policy(policy.PolicyString(
            """#<?cfg paf policy?>
            SHAPELET_MODEL.enabled: true
            """ 
            ))
        self.mp.configure(self.pol)        

    def tearDown(self):
        del self.exp
        del self.source
        del self.mp
        del self.pol

    def testBadInput(self):
        #no Source
        meas = self.mp.measure(afwDetection.Peak(), None)        
        photom = meas.find("SHAPELET_MODEL")
        status = int(photom.get("status"))
        self.assertTrue(status & measAlgorithms.Flags.PHOTOM_NO_SOURCE)

        #bad sources
        noFp = afwDetection.Source()
        noFp.setXAstrom(self.source.getXAstrom())
        noFp.setYAstrom(self.source.getYAstrom())
        noFp.setIxx(self.source.getIxx())
        noFp.setIyy(self.source.getIyy())
        noFp.setIxy(self.source.getIxy())
        meas=self.mp.measure(afwDetection.Peak(), noFp)
        photom = meas.find("SHAPELET_MODEL")
        status = int(photom.get("status"))
        self.assertTrue(status & measAlgorithms.Flags.PHOTOM_NO_FOOTPRINT)

        badMoments = afwDetection.Source(self.source)
        badMoments.setIxx(float('nan'))
        meas=self.mp.measure(afwDetection.Peak(), badMoments)
        photom = meas.find("SHAPELET_MODEL")
        status = int(photom.get("status"))
        self.assertTrue(status & measAlgorithms.Flags.SHAPELET_PHOTOM_BAD_MOMENTS)

        #badMoments.setIxx(float('inf'))
        #meas=self.mp.measure(afwDetection.Peak(), badMoments)
        #photom = meas.find("SHAPELET_MODEL")
        #status = int(photom.get(afwDetection.Schema("STATUS", 2, afwDetection.Schema.INT)))
        #self.assertTrue(status & 0x020)

        badMoments.setIxx(100000.0)
        meas=self.mp.measure(afwDetection.Peak(), badMoments)
        photom = meas.find("SHAPELET_MODEL")
        status = int(photom.get("status"))
        self.assertTrue(status & measAlgorithms.Flags.SHAPELET_PHOTOM_BAD_MOMENTS)

        #no psf
        badExp = self.exp.Factory(self.exp, True)
        badExp.setPsf(None)

        mp = measAlgorithms.makeMeasurePhotometry(badExp)
        mp.configure(self.pol)
        mp.addAlgorithm("SHAPELET_MODEL")

        meas =mp.measure(afwDetection.Peak(), self.source)        
        photom = meas.find("SHAPELET_MODEL")
        status = int(photom.get("status"))
        self.assertTrue(status & measAlgorithms.Flags.PHOTOM_NO_PSF)
   
        #No pixels
        badExp.setPsf(self.exp.getPsf())
        mask = badExp.getMaskedImage().getMask()
        mask |= afwImage.MaskU.getPlaneBitMask("BAD")
        
        mp.setImage(badExp)
        meas =mp.measure(afwDetection.Peak(), self.source)        
        photom = meas.find("SHAPELET_MODEL")
        status = int(photom.get("status"))
        self.assertTrue(status & measAlgorithms.Flags.PHOTOM_NO_FOOTPRINT)

    def testPsfFlux(self):
        """Test that fluxes are measured correctly"""
        #
        # Various algorithms
        #
        meas = self.mp.measure(afwDetection.Peak(self.xc, self.yc), self.source)
        print "SHAPELET_MODEL"
        photom = meas.find("SHAPELET_MODEL")
        print >> sys.stderr, "flux:", photom.get("flux")
        print >> sys.stderr, "fluxErr:",photom.get("fluxErr")
        print >> sys.stderr, "status:", int(photom.get("status"))
        print >> sys.stderr, "e1:",photom.get("e1")
        print >> sys.stderr, "e2:", photom.get("e2")
        print >> sys.stderr, "radius",photom.get("radius")

class MultifitDataTestCase(unittest.TestCase):
    def setUp(self):
        bf = ButlerFactory(mapper=lsst.meas.multifitData.DatasetMapper())
        self.butler = bf.create()
        self.pol = policy.Policy(policy.PolicyString(
            """#<?cfg paf policy?>
            SHAPELET_MODEL.enabled: true
            """ 
            ))

    def tearDown(self):
        del self.butler
        del self.pol

    def testAll(self):
        for i in range(10):
            exp = self.butler.get('exp', id=i)
            psf = self.butler.get('psf', id=i)
            exp.setPsf(psf)
            sources = self.butler.get('src', id=i)

            mp = measAlgorithms.makeMeasurePhotometry(exp)
            mp.addAlgorithm("SHAPELET_MODEL")
            mp.configure(self.pol)

            for s in sources:
                mp.measure(afwDetection.Peak(), s)




#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilTests.init()
    
    suites = []
    suites += unittest.makeSuite(GaussianPsfTestCase)
    #if (eups.productDir("meas_multifitData")):
    #    suites += unittest.makeSuite(MultifitDataTestCase)
    #else:
    #    print "meas_multifitData is not setup, skipping extensive tests"
    suites += unittest.makeSuite(utilTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
