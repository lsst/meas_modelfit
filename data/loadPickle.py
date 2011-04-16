import cPickle
import lsst.meas.multifit

s = cPickle.load(open("sersic-basis.p"))

for k, v in s.iteritems():
    items, forward, reverse = v
    components = lsst.meas.multifit.CompoundShapeletBuilder.ComponentVector(
        [lsst.meas.multifit.ShapeletModelBasis.make(order, scale) for order, scale in items]
        )
    builder = lsst.meas.multifit.CompoundShapeletBuilder(components)
    builder.setMapping(forward, reverse)
    basis = builder.build()
    basis.save("%s.boost" % k)
