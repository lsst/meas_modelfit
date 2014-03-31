import lsst.meas.multifit
import lsst.meas.extensions.multiShapelet

root.measurement.algorithms.names |= ["multishapelet.psf", "cmodel"]
root.measurement.slots.modelFlux = "cmodel.flux"
