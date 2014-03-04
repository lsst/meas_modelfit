import os
import lsst.daf.persistence
from matplotlib import pyplot

FLUXES = [("flux.psf", "mag.psf"), ("flux.kron", "mag.kron"), ("flux.sinc", "mag.sinc"),
          ("cmodel.flux", "cmodel.mag"), ("cmodel.initial.flux", "cmodel.initial.mag"),
          ("cmodel.exp.flux", "cmodel.exp.mag"), ("cmodel.dev.flux", "cmodel.dev.mag"),
          ]

def buildCatalog(self, root=None, rerun=None, butler=None, noCache=False,
                 doPrintStatus=False, doPrintMissing=False, prefix="",
                 **dataId):
    if butler is None:
        if root is None:
            root = os.path.join(os.environ("SUPRIME_DATA_DIR"), "rerun", rerun)
        butler = lsst.daf.persistence.Butler(root)
    schemaDataSet = prefix + "src_schema"
    inputSchema = butler.get(schemaDataSet, immediate=True).schema
    mapper = lsst.afw.table.SchemaMapper(inputSchema)
    mapper.addMinimalSchema(inputSchema)
    visitKey = mapper.addOutputField(lsst.afw.table.Field["I"]("visit", "visit number of exposure"))
    ccdKey = mapper.addOutputField(lsst.afw.table.Field["I"]("ccd", "ccd number of exposure"))
    magKeys = []
    for inputName, outputName in FLUXES:
        magKey = mapper.addOutputField(
            lsst.afw.Field["D"](outputName, "magnitude for %s" % inputName, "mag")
            )
        magErrKey = mapper.addOutputField(
            lsst.afw.Field["D"](outputName + ".err", "magnitude error for %s" % inputName, "mag")
            )
        magKeys.append((inputSchema.find(inputName).key, inputSchema.find(inputName + ".err").key,
                        magKey, magErrKey))
    srcDataSet = prefix + "src"
    calexpDataSet = prefix + "calexp"
    srcCache = {}
    fullSize = 0
    for dataRef in butler.subset(calexpDataSet, **dataId):
        if not dataRef.datasetExists(srcDataSet):
            if doPrintMissing:
                print "%s missing for data id %s" % (srcDataSet, dataRef.dataId)
            continue

        if doPrintStatus:
            print "Loading catalog for %s" % dataRef.dataId

        src = dataRef.get(srcDataSet, immediate=True)
        if not noCache:
            srcCache[dataRef.dataId] = src

        fullSize += len(src)

    outCat = lsst.afw.table.SourceCatalog(mapper.getOutputSchema())
    outCat.reserve(fullSize)
    for dataRef in butler.subset(calexpDataSet, **dataId):
        if not dataRef.datasetExists(srcDataSet):
            continue

        if doPrintStatus:
            print "Calibrating and appending catalog for %s" % dataRef.dataId

        if noCache:
            src = dataRef.get(srcDataSet, immediate=True)
        else:
            src = srcCache.pop(dataRef.dataId)

        # transfer the all src fields to the output catalog, and create a view
        # into the output catalog containing only the new records
        n1 = len(outCat)
        outCat.extend(src, mappper=mapper)
        del src
        n2 = len(outCat)
        subOutCat = outCat[n1:n2]

        # Set the data ID fields in the new records
        subOutCat[ccdKey] = dataId.ccd
        subOutCat[visitKey] = dataId.visit

        # Set the magnitude fields in the new records
        calexp = dataRef.get(calexpDataSet, immediate=True)
        calib = calexp.getCalib()
        del calexp
        calib.setThrowOnNegativeFlux(False)
        for fluxKey, fluxErrKey, magKey, magErrKey in magKeys:
            subOutCat[magKey], subOutCat[magErrKey] = calib.getMagnitude(subOutCat[fluxKey],
                                                                         subOutCat[fluxErrKey])
    return outCat

def getFlagMask(catalog, *flags):
    return numpy.logical_not(numpy.logical_or.reduce([catalog[flag] for flag in flags]))

def plotDelta(catalog, x, y, c=None, **kwds):
    kwds.setdefault("alpha", 0.2)
    kwds.setdefault("linewidth", 0.0)
    kwds.setdefault("s", 4)
    pyplot.scatter(catalog[x], catalog[y] - catalog[x], c=(None if c is None else catalog[c]), **kwds)
