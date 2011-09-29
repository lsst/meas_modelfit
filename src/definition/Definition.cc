#include "lsst/meas/multifit/definition/Definition.h"

namespace lsst { namespace meas { namespace multifit { namespace definition {

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), _wcs(other._wcs) {}

}}}} // namespace lsst::meas::multifit::definition
