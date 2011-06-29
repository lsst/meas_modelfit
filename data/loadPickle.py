import cPickle
import lsst.meas.multifit
import numpy

s = cPickle.load(open("sersic-basis.p"))

for k, v in s.iteritems():
    items, forward, reverse = v
    components = lsst.meas.multifit.CompoundShapeletBuilder.ComponentVector(
        [lsst.meas.multifit.ShapeletModelBasis.make(order, scale) for order, scale in items]
        )
    builder = lsst.meas.multifit.CompoundShapeletBuilder(components)
    builder.setMapping(forward, reverse)
    builder.normalize()
    #constraint = numpy.zeros((min(3, builder.getSize()), builder.getSize()), dtype=float)
    #constraint[:,:constraint.shape[0]] = numpy.identity(constraint.shape[0], dtype=float)
    #builder.setConstraint(constraint, numpy.zeros(constraint.shape[0], dtype=float))
    basis = builder.build()
    basis.save("%s.boost" % k)
