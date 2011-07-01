#include "lsst/meas/multifit/grid/FluxGroup.h"
#include <boost/format.hpp>

namespace lsst { namespace meas { namespace multifit { namespace grid {

void FluxGroup::initialize() {
    int constraintCount = 0;
    for (ComponentArray::iterator i = components.begin(); i != components.end(); ++i) {
        if (i->getBasis()) {
            constraintCount += i->getBasis()->getConstraintSize();
        } else {
            constraintCount += 1;
        }
    }
    _constraintMatrix = ndarray::allocate(getSourceCoefficientCount(), constraintCount);
    _constraintMatrix.deep() = 0.0;
    _constraintVector = ndarray::allocate(constraintCount);
    _constraintVector.deep() = 0.0;
    _integration = ndarray::allocate(getSourceCoefficientCount());
    _integration.deep() = 0.0;
    int coefficientOffset = 0;
    int constraintOffset = 0;
    for (ComponentArray::iterator i = components.begin(); i != components.end(); ++i) {
        if (i->getBasis()) {
            _constraintMatrix[
                ndarray::view(
                    coefficientOffset, coefficientOffset + i->getBasis()->getSize()
                )(
                    constraintOffset, constraintOffset + i->getBasis()->getConstraintSize()
                )
            ].deep() = i->getBasis()->getConstraintMatrix();
            _constraintVector[
                ndarray::view(constraintOffset, constraintOffset + i->getBasis()->getConstraintSize())
            ].deep() = i->getBasis()->getConstraintVector();
            constraintOffset += i->getBasis()->getConstraintSize();
            i->getBasis()->integrate(
                _integration[ndarray::view(coefficientOffset, coefficientOffset + i->getBasis()->getSize())]
            );
        } else {
            _constraintMatrix[coefficientOffset][constraintOffset] = 1.0;
            // constraint vector for point sources is 0
            constraintOffset += 1;
            _integration[coefficientOffset] = 1.0;
        }
    }
}

}}}} // namespace lsst::meas::multifit::grid
