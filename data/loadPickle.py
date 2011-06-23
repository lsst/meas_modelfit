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
    builder.setMapping(forward)
    if builder.getSize() > 2:
        assert(builder.getSize() != 3)
        constraint = numpy.zeros((4, builder.getSize()), dtype=float)
        constraint[:3,:3] = numpy.identity(3, dtype=float)
        integral = numpy.zeros(builder.getSize(), dtype=float)
        builder.integrate(integral)
        constraint[3, 3:] = integral[3:] / (integral[3:]**2).sum()**0.5
    else:
        constraint = numpy.identity(2, dtype=float)
    builder.setConstraint(constraint, numpy.zeros(constraint.shape[0], dtype=float))
    basis = builder.build()
    basis.save("%s.boost" % k)
