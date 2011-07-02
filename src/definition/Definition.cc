#include "lsst/meas/multifit/definition/Definition.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/ndarray/eigen.h"
#include <limits>
#include <Eigen/Array>

namespace lsst { namespace meas { namespace multifit { namespace definition {

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), _wcs(other._wcs) {}

}}}} // namespace lsst::meas::multifit::definition
