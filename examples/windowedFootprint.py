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
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult
import numpy


def main():
    fp = afwDet.Footprint(afwImg.BBox(afwImg.PointI(0,0), 2, 3))
    fp.normalize()
    window = afwGeom.BoxI(afwGeom.makePointI(1,0), afwGeom.makeExtentI(1,3))
    wfp = measMult.WindowedFootprint(fp, window)
    full = numpy.array([[1., 2.], [3., 4.], [5., 6.]])
    compressed = numpy.zeros(fp.getNpix())
    wfp.compress(full, compressed)

    print full
    print window
    print compressed

if __name__ == '__main__':
    main()
    

