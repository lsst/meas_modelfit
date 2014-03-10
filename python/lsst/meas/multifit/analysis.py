import os
import numpy
from matplotlib import pyplot

import lsst.obs.hsc
import lsst.obs.suprimecam
import lsst.daf.persistence

FLUXES = ["flux.psf", "flux.kron", "flux.sinc",
          "cmodel.flux", "cmodel.initial.flux", "cmodel.exp.flux", "cmodel.dev.flux",
          ]

def getRawDataIds(*visits):
    r = []
    for visit in visits:
        r.extend(dict(ccd=ccd, visit=visit) for ccd in xrange(0, 104))
    return r

def getCoaddDataIds(filter, patchesX, patchesY, tract=0):
    r = []
    for patchX in patchesX:
        for patchY in patchesY:
            r.append(dict(filter=filter, tract=tract, patch="%d,%d" % (patchX, patchY)))
    return r

class CalculatedFieldSet(object):
    """A base class for objects that allocate and fill one or more calculated catalog columns."""

    def fillCatalog(self, catalog, dataRef, calexp):
        pass

    def fillRecord(self, record, dataRef, calexp):
        pass

class CModelLogRadius(CalculatedFieldSet):
    """A CalculatedFieldSet that computes the log10 of the effective radius in arcseconds
    for the CModel algorithm, using the fracDev to weigh the exponential and de Vaucouleur
    radii."""

    def __init__(self, schema):
        self.outKey = schema.addField("cmodel.log10r", doc="logarithm of half-light radius",
                                      units="log(arcsec)", type=float)
        self.rExpKey = schema.find("cmodel.exp.nonlinear").key[2]
        self.rDevKey = schema.find("cmodel.dev.nonlinear").key[2]
        self.fracDevKey = schema.find("cmodel.fracDev").key
        self.logFactor = 1.0 / numpy.log(10.0)

    def fillCatalog(self, catalog, dataRef, calexp):
        catalog[self.outKey][:] = (catalog[self.rExpKey] * (1.0 - catalog[self.fracDevKey])
                                   + catalog[self.rDevKey] * catalog[self.fracDevKey]) * self.logFactor

    def fillRecord(self, record, dataRef, calexp):
        pass

class CModelEllipticity(CalculatedFieldSet):
    """A CalculatedFieldSet that computes the log10 of the effective radius in arcseconds
    for the CModel algorithm, using the fracDev to weigh the exponential and de Vaucouleur
    radii."""

    def __init__(self, schema):
        self.outKey = schema.addField("cmodel.ellipticity", doc="logarithm of half-light radius",
                                      units="log(arcsec)", type=float)
        self.e1ExpKey = schema.find("cmodel.exp.nonlinear").key[0]
        self.e1DevKey = schema.find("cmodel.dev.nonlinear").key[0]
        self.e2ExpKey = schema.find("cmodel.exp.nonlinear").key[1]
        self.e2DevKey = schema.find("cmodel.dev.nonlinear").key[1]
        self.fracDevKey = schema.find("cmodel.fracDev").key

    def fillCatalog(self, catalog, dataRef, calexp):
        eExp = (catalog[self.e1ExpKey]**2 + catalog[self.e2ExpKey]**2)**0.5
        eDev = (catalog[self.e1DevKey]**2 + catalog[self.e2DevKey]**2)**0.5
        catalog[self.outKey][:] = (eExp * (1.0 - catalog[self.fracDevKey])
                                   + eDev * catalog[self.fracDevKey])

    def fillRecord(self, record, dataRef, calexp):
        pass

class Magnitudes(CalculatedFieldSet):
    """A CalculatedFieldSet that computes magnitudes corresponding to all flux fields
    in the module-scope variable FLUXES."""

    def __init__(self, schema):
        self.keyDicts = []
        for flux in FLUXES:
            mag = flux.replace("flux", "mag")
            keys = {"flux": schema.find(flux).key, "flux.err":schema.find(flux+".err").key}
            keys["mag"] = schema.addField(mag, doc="magnitude for %s" % flux, units="mag", type=float)
            keys["mag.err"] = schema.addField(mag + ".err", doc="magnitude error for %s" % flux,
                                              units="mag", type=float)
            self.keyDicts.append(keys)

    def fillCatalog(self, catalog, dataRef, calexp):
        calib = calexp.getCalib()
        calib.setThrowOnNegativeFlux(False)
        for keys in self.keyDicts:
            catalog[keys["mag"]][:], catalog[keys["mag.err"]][:] = \
                calib.getMagnitude(catalog[keys["flux"]], catalog[keys["flux.err"]])

FORCED_FILTERS = {'HSC-G': 'g', 'HSC-I': 'i', 'HSC-Z': 'z'}

class ForcedMagnitudes(object):
    """An object that computes magnitudes from associated forced catalogs
    for all fluxes in the module-scope variable FLUXES."""

    def __init__(self, outputSchema, forcedSchema):
        self.keys = {}
        self.objectIdKey = forcedSchema.find("object.id").key
        for dataIdFilter, fieldNameFilter in FORCED_FILTERS.iteritems():
            prefix = fieldNameFilter + "."
            self.keys[dataIdFilter] = []
            for flux in FLUXES:
                mag = flux.replace("flux", "mag")
                algKeys = {"flux": forcedSchema.find(flux).key,
                             "flux.err": forcedSchema.find(flux + ".err").key}
                algKeys["mag"] = outputSchema.addField(prefix + mag,
                                                       doc="magnitude for %s" % flux,
                                                       units="mag", type=float)
                algKeys["mag.err"] = outputSchema.addField(prefix + mag + ".err",
                                                           doc="magnitude error for %s" % flux,
                                                           units="mag", type=float)
                algKeys["flag"] = outputSchema.addField(prefix + mag + ".flags", type="Flag",
                                                       doc="flag for forced photometry")
                self.keys[dataIdFilter].append(algKeys)

    def fillCatalog(self, catalog, dataRef, calexp):
        for dataIdFilter, filterKeys in self.keys.iteritems():
            dataId = dataRef.dataId.copy()
            dataId["filter"] = dataIdFilter
            forcedDataRef = dataRef.butlerSubset.butler.dataRef("deepCoadd_forced_src", **dataId)
            try:
                forcedSrc = forcedDataRef.get("deepCoadd_forced_src", immediate=True)
            except Exception:
                for record in catalog:
                    for algKeys in filterKeys:
                        record.set(algKeys["flag"], True)
                continue
            forcedCalib = forcedDataRef.get("deepCoadd_calexp", immediate=True).getCalib()
            forcedCalib.setThrowOnNegativeFlux(False)
            forcedSrc.sort(self.objectIdKey)
            forcedSrcIter = iter(forcedSrc)
            for record in catalog:
                try:
                    forcedSrcRecord = forcedSrcIter.next()
                except StopIteration:
                    pass
                if forcedSrcRecord.get(self.objectIdKey) == record.getId():
                    for algKeys in filterKeys:
                        mag, magErr = forcedCalib.getMagnitude(forcedSrcRecord.get(algKeys["flux"]),
                                                               forcedSrcRecord.get(algKeys["flux.err"]))
                        record.set(algKeys["mag"], mag)
                        record.set(algKeys["mag.err"], magErr)
                elif forcedSrcRecord.get(self.objectIdKey) > record.getId():
                    for algKeys in filterKeys:
                        record.set(algKeys["flag"], True)
                else:
                    while forcedSrcRecord.get(self.objectIdKey) < record.getId():
                        try:
                            forcedSrcRecord = forcedSrcIter.next()
                        except StopIteration:
                            pass

class SubaruDataIdSchemaManager(CalculatedFieldSet):
    """A CalculatedFieldSet appropriate for managing Data ID fields for non-coadd Subaru data"""

    def __init__(self, schema):
        self.visitKey = schema.addField(lsst.afw.table.Field["I"]("visit", "visit number of exposure"))
        self.ccdKey = schema.addField(lsst.afw.table.Field["I"]("ccd", "ccd number of exposure"))
        self.filterKey = schema.addField(lsst.afw.table.Field[str]("filter", "filter name", 8))

    def fillRecord(self, record, dataRef, calexp):
        record.set(self.filterKey, dataRef.dataId["filter"])
        record.set(self.ccdKey, dataRef.dataId["ccd"])
        record.set(self.visitKey, dataRef.dataId["visit"])

class CoaddDataIdSchemaManager(CalculatedFieldSet):
    """A CalculatedFieldSet appropriate for managing Data ID fields for coadd data"""

    def __init__(self, schema):
        self.tractKey = schema.addField(lsst.afw.table.Field["I"]("tract", "coadd tract"))
        self.patchKey = schema.addField(lsst.afw.table.Field[str]("patch", "coadd patch", 32))
        self.filterKey = schema.addField(lsst.afw.table.Field[str]("filter", "filter name", 32))

    def fillCatalog(self, catalog, dataRef, calexp):
        pass

    def fillRecord(self, record, dataRef, calexp):
        record.set(self.patchKey, dataRef.dataId["patch"])
        record.set(self.tractKey, dataRef.dataId["tract"])
        record.set(self.filterKey, dataRef.dataId["filter"])


def makeDataIdSchemaManager(schema, butler, prefix):
    """Return a CalculatedFieldSet that manages Data ID fields appropriate for the given
    butler and data prefix."""
    if prefix.endswith("Coadd_"):
        return CoaddDataIdSchemaManager(schema)
    elif prefix == "" and (
        isinstance(butler.mapper, lsst.obs.hsc.HscMapper)
        or isinstance(butler.mapper, lsst.obs.suprimecam.SuprimecamMapperBase)
        ):
        return SubaruDataIdSchemaManager(schema)
    else:
        raise ValueError("Unknown mapper")

CALCULATED = (CModelLogRadius, CModelEllipticity)

def printMissing(dataSet, root=None, rerun=None, butler=None, dataIds=[]):
    if butler is None:
        if root is None:
            root = os.path.join(os.environ["SUPRIME_DATA_DIR"], "rerun", rerun)
        butler = lsst.daf.persistence.Butler(root)
    for dataId in dataIds:
        if not butler.datasetExists(dataSet, dataId):
            print dataId

def buildCatalog(
    root=None, rerun=None, butler=None, noCache=False, dataIds=[],
    doPrintStatus=False, doPrintMissing=False, prefix="",
    doReadForcedCoadd=False, calculated=CALCULATED
    ):

    if butler is None:
        if root is None:
            root = os.path.join(os.environ["SUPRIME_DATA_DIR"], "rerun", rerun)
        butler = lsst.daf.persistence.Butler(root)
    schemaDataSet = prefix + "src_schema"
    inputSchema = butler.get(schemaDataSet, immediate=True).schema
    mapper = lsst.afw.table.SchemaMapper(inputSchema)
    mapper.addMinimalSchema(inputSchema)

    # Setup calculated fields, including Data ID fields
    calculated = [cls(mapper.editOutputSchema()) for cls in calculated]
    calculated.append(Magnitudes(mapper.editOutputSchema()))
    calculated.append(makeDataIdSchemaManager(mapper.editOutputSchema(), butler, prefix))

    if doReadForcedCoadd:
        forcedDataSet = prefix + "forced_src"
        outputSchema = mapper.editOutputSchema()
        forcedSchema = butler.get(prefix + "forced_src_schema", immediate=True).schema
        calculated.append(ForcedMagnitudes(outputSchema, forcedSchema))

    srcDataSet = prefix + "src"
    calexpDataSet = prefix + "calexp"
    srcCache = {}
    fullSize = 0

    for dataId in dataIds:
        dataRef = butler.dataRef(srcDataSet, **dataId)
        dataId = dataRef.dataId

        if not dataRef.datasetExists(srcDataSet):
            if doPrintMissing:
                print "%s missing for data id %s" % (srcDataSet, dataId)
            continue

        if doPrintStatus:
            print "Loading catalog for %s" % dataId

        src = dataRef.get(srcDataSet, flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS, immediate=True)
        if not noCache:
            srcCache[tuple(dataId.itervalues())] = src

        fullSize += len(src)

    outCat = lsst.afw.table.SourceCatalog(mapper.getOutputSchema())
    outCat.reserve(fullSize)
    for dataId in dataIds:
        dataRef = butler.dataRef(srcDataSet, **dataId)
        dataId = dataRef.dataId
        if noCache:
            if not dataRef.datasetExists(srcDataSet):
                continue
            src = dataRef.get(srcDataSet, flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS, immediate=True)
        else:
            src = srcCache.pop(tuple(dataId.itervalues()), None)
            if src is None:
                continue

        if doPrintStatus:
            print "Calibrating and appending catalog for %s" % dataId

        # transfer the all src fields to the output catalog, and create a view
        # into the output catalog containing only the new records
        n1 = len(outCat)
        outCat.extend(src, mapper=mapper)
        del src
        n2 = len(outCat)
        subOutCat = outCat[n1:n2]

        # Fill custom calculated fields
        calexp = butler.get(calexpDataSet, dataId, immediate=True)
        for field in calculated:
            field.fillCatalog(subOutCat, dataRef, calexp)
        del calexp

    return outCat

def getFlagMask(catalog, *flags, **kwds):
    bad = [catalog[flag] for flag in flags]
    if kwds.get("isolated", False):
        bad.append(catalog["parent"] != 0)
        bad.append(catalog["deblend.nchild"] != 0)
    elif kwds.get("deblended", False):
        bad.append(catalog["deblend.nchild"] != 0)
    return numpy.logical_not(numpy.logical_or.reduce(bad))

def plotDelta(catalog, x, y, c=None, **kwds):
    kwds.setdefault("alpha", 0.2)
    kwds.setdefault("linewidth", 0.0)
    kwds.setdefault("s", 4)
    pyplot.scatter(catalog[x], catalog[y] - catalog[x], c=(None if c is None else catalog[c]), **kwds)

def plotScatter(catalog, x, y, c=None, **kwds):
    kwds.setdefault("alpha", 0.2)
    kwds.setdefault("linewidth", 0.0)
    kwds.setdefault("s", 4)
    pyplot.scatter(catalog[x], catalog[y], c=(None if c is None else catalog[c]), **kwds)
