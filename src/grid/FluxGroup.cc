#include "lsst/meas/multifit/grid/FluxGroup.h"
#include <boost/format.hpp>

namespace lsst { namespace meas { namespace multifit { namespace grid {

void FluxGroup::initialize() {
    int constraintCount = 0;
    for (ComponentArray::iterator i = components.begin(); i != components.end(); ++i) {
        if (i->getBasis()) {
            constraintCount += i->getBasis()->getConstraintCount();
        } else {
            constraintCount += 1;
        }
    }
    _constraintMatrix = ndarray::allocate(constraintCount, getSourceCoefficientCount());
    _constraintMatrix.deep() = 0.0;
    _constraintVector = ndarray::allocate(constraintCount);
    _constraintVector.deep() = 0.0;
    _integration = ndarray::allocate(getSourceCoefficientCount());
    _integration.deep() = 0.0;
    int coefficientOffset = 0;
    int constraintOffset = 0;
    for (ComponentArray::iterator i = components.begin(); i != components.end(); ++i) {
        if (i->getBasis()) {
            if (i->getBasis()->getConstraintCount() > 0) {
                _constraintMatrix[
                    ndarray::view(
                        constraintOffset, constraintOffset + i->getBasis()->getConstraintCount()
                    )(
                        coefficientOffset, coefficientOffset + i->getBasis()->getSize()
                    )
                ] = i->getBasis()->getConstraintMatrix();
                _constraintVector[
                    ndarray::view(constraintOffset, constraintOffset + i->getBasis()->getConstraintCount())
                ] = i->getBasis()->getConstraintVector();
            }
            _integration[ndarray::view(coefficientOffset, coefficientOffset + i->getBasis()->getSize())]
                = i->getBasis()->getIntegration();
            constraintOffset += i->getBasis()->getConstraintCount();
            coefficientOffset += i->getBasis()->getSize();
        } else {
            _constraintMatrix[constraintOffset][coefficientOffset] = 1.0;
            // constraint vector for point sources is 0
            _integration[coefficientOffset] = 1.0;
            constraintOffset += 1;
            coefficientOffset += 1;
        }
    }
}

}}}} // namespace lsst::meas::multifit::grid
