#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

int const Grid::getFilterIndex(FilterId filterId) const {
    FilterMap::const_iterator i = _filters.find(filterId);
    if(i == _filters.end()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Filter with ID %d not found.") % filterId).str()
        );
    }
    return i->second;
}

void Grid::writeParameters(ndarray::Array<double,1,1> const & parameters) const {
    double * paramIter = parameters.getData();
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        detail::SharedElementTraits<POSITION>::writeIter(paramIter + i->offset, i->getValue());
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        detail::SharedElementTraits<RADIUS>::writeIter(paramIter + i->offset, i->getValue());
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        detail::SharedElementTraits<ELLIPTICITY>::writeIter(paramIter + i->offset, i->getValue());
    }
}

bool Grid::checkBounds(ndarray::Array<double const,1,1> const & parameters) const {
    double const * paramIter = parameters.getData();
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        if (!i->getBounds().checkBounds(paramIter + i->offset)) return false;
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        if (!i->getBounds().checkBounds(paramIter + i->offset)) return false;
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        if (!i->getBounds().checkBounds(paramIter) + i->offset) return false;
    }
    return true;
}

double Grid::clipToBounds(ndarray::Array<double,1,1> const & parameters) const {
    double * paramIter = parameters.getData();
    double value = 0.0;
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        value += i->getBounds().clipToBounds(paramIter + i->offset);
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        value += i->getBounds().clipToBounds(paramIter + i->offset);
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        value += i->getBounds().clipToBounds(paramIter + i->offset);
    }
    return value;
}

}}}} // namespace lsst::meas::multifit::grid
