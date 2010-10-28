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

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.multifit as measMult
import numpy
import numpy.random

from makeImageStack import makeImageStack

def initializeModelEvaluator():
    flux = 1.0
    centroid = afwGeom.makePointD(45, 45)
    psModel = measMult.createPointSourceModel(flux, centroid)

    (miList, psfList, transformList) = makeImageStack(psModel, 5)
    modelEvaluator = measMult.ModelEvaluator(psModel)
    modelEvaluator.setData(miList, psfList, transformList)

    numpy.set_printoptions(threshold=numpy.nan)
    print "ModelEvaluator nProjections: %d"%modelEvaluator.getNProjections()
    print "ModelEvaluator nPixels: %d"%modelEvaluator.getNPixels()

    print "ModelEvaluator image vector: %s"%modelEvaluator.getDataVector()
    print "ModelEvaluator variance vector: %s"%modelEvaluator.getSigmaVector()  
    print "ModelEvaluator ModelImage: %s"%modelEvaluator.computeModelImage()
    print "ModelEvaluator LinearParameterDerivative: %s"%modelEvaluator.computeLinearParameterDerivative()
    print "ModelEvaluator NoninearParameterDerivative: %s"%modelEvaluator.computeNonlinearParameterDerivative()


if __name__ == "__main__":
    initializeModelEvaluator()
