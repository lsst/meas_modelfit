import lsst.meas.multifitData 
import lsst.afw.detection
import lsst.afw.image
import lsst.meas.multifit
import lsst.meas.algorithms
import lsst.daf.persistence
import numpy
import sys
import logging

try:
    import cPickle as pickle
except ImportError:
    import pickle

FLUX_SCHEMA = lsst.afw.detection.Schema("Flux", 0, lsst.afw.detection.Schema.DOUBLE)
FLUX_ERR_SCHEMA = lsst.afw.detection.Schema("FluxErr", 1, lsst.afw.detection.Schema.DOUBLE)
STATUS_SCHEMA = lsst.afw.detection.Schema("Status", 2, lsst.afw.detection.Schema.INT)
E1_SCHEMA = lsst.afw.detection.Schema("E1", 3, lsst.afw.detection.Schema.DOUBLE)
E2_SCHEMA = lsst.afw.detection.Schema("E2", 4, lsst.afw.detection.Schema.DOUBLE)
RADIUS_SCHEMA = lsst.afw.detection.Schema("R", 5, lsst.afw.detection.Schema.DOUBLE)
COEFF_SCHEMA = lsst.afw.detection.Schema("COEFFICIENTS", 6, lsst.afw.detection.Schema.DOUBLE)

fields = (("dataset", int),
          ("id", numpy.int64),
          ("psf_flux", float), 
          ("psf_flux_err", float),
          ("x", float),
          ("y", float),
          ("ixx", float),
          ("iyy", float),
          ("ixy", float),
          ("status", numpy.int64),
          ("flux", float),
          ("flux_err", float),
          ("e1", float),
          ("e2", float),
          ("r", float),
          ("pixels", int),
          )


def fit(datasets=(0,1,2,3,4,5,6,7,8,9), basisSize=8, fitDelta=True):
    bf = lsst.daf.persistence.ButlerFactory( \
        mapper=lsst.meas.multifitData.DatasetMapper())
    butler = bf.create()
    algorithm = "SHAPELET_MODEL_%d" % basisSize
    policy = lsst.pex.policy.Policy()
    policy.add(algorithm + ".enabled", True)
    policy.add(algorithm + ".fitDeltaFunction", fitDelta)
    nCoeffs = basisSize
    if fitDelta:
        nCoeffs += 1
    if basisSize == 2:
        basis = lsst.meas.multifit.utils.loadBasis("ed+02:0000")
    elif basisSize == 8:
        basis = lsst.meas.multifit.utils.loadBasis("ed+06:2000")
    elif basisSize == 17:
        basis = lsst.meas.multifit.utils.loadBasis("ed+15:4000")
    rawIntegral = numpy.zeros(nCoeffs, dtype=float)
    basis.integrate(rawIntegral[:basisSize])
    if fitDelta:
        rawIntegral[-1] = 1.0
    dtype = numpy.dtype(list(fields) + [("integral", float, nCoeffs), ("coeff", float, nCoeffs)])
    tables = []
    for d in datasets:
        logging.info("Processing dataset %d" % d)
        if not butler.datasetExists("src", id=d):
            continue
        psf = butler.get("psf", id=d)
        exp = butler.get("exp", id=d)
        exp.setPsf(psf)
        sources = butler.get("src", id=d)

        measurePhotometry = lsst.meas.algorithms.makeMeasurePhotometry(exp)
        measurePhotometry.addAlgorithm(algorithm)
        measurePhotometry.configure(policy)

        table = numpy.zeros((len(sources)), dtype=dtype)
        table["dataset"] = d
        for j, src in enumerate(sources):
            logging.debug("Fitting source %d (%d/%d)" % (src.getId(), j+1, len(sources)))
            record = table[j]
            record["id"] = src.getId()
            record["psf_flux"] = src.getPsfFlux()
            record["psf_flux_err"] = src.getPsfFluxErr()
            record["x"] = src.getXAstrom()
            record["y"] = src.getYAstrom()
            record["ixx"] = src.getIxx()
            record["iyy"] = src.getIyy()
            record["ixy"] = src.getIxy()
            record["pixels"] = src.getFootprint().getArea()
            photom = measurePhotometry.measure(lsst.afw.detection.Peak(), src).find(algorithm)
            record["status"] = photom.get(STATUS_SCHEMA)
            record["flux"] = photom.get(FLUX_SCHEMA)
            record["flux_err"] = photom.get(FLUX_ERR_SCHEMA)
            record["e1"] = photom.get(E1_SCHEMA)
            record["e2"] = photom.get(E2_SCHEMA)
            record["r"] = photom.get(RADIUS_SCHEMA)
            ellipse = lsst.meas.multifit.EllipseCore(record["e1"], record["e2"], record["r"])
            f = ellipse.getArea() / numpy.pi
            record["integral"][:] = rawIntegral
            record["integral"][:basisSize] *= f
            for i in range(basisSize):
                record["coeff"][i] = photom.get(i, COEFF_SCHEMA)
            if fitDelta:
                record["coeff"][basisSize] = photom.get(basisSize, COEFF_SCHEMA)
            altFlux = numpy.dot(record["integral"], record["coeff"])
            if not numpy.allclose(altFlux, record["flux"]):
                logging.warning("altFlux=%s, flux=%s" % (altFlux, record["flux"]))
        tables.append(table)

    return numpy.concatenate(tables)

def build(filename, *args, **kwds):
    logging.basicConfig(level=logging.DEBUG)
    full = fit(*args, **kwds)
    outfile = open(filename, 'wb')
    pickle.dump(full, outfile, protocol=2)
    outfile.close()

def load(filename, filter_status=True, filter_flux=True, dataset=None):
    infile = open(filename, "rb")
    table = pickle.load(infile)
    if filter_status:
        table = table[table["status"] == 0]
    if filter_flux:
        table = table[table["flux"] > 0]
    if dataset is not None:
        table = table[table["dataset"] == dataset]
    return table

def m_radius(table):
    return (table["ixx"] + table["iyy"])**0.5

def ellipticity(table):
    return (table["e1"]**2 + table["e2"]**2)**0.5

