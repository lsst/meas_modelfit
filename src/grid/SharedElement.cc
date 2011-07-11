#include "lsst/meas/multifit/grid/SharedElement.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

template <SharedElementType E>
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & component) {
    return os << (*definition::SharedElement<E>::make(component.getValue(), component.isActive()))
              << " [" << component.offset << "]";
}

template std::ostream & operator<<(std::ostream &, SharedElement<POSITION> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<RADIUS> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<ELLIPTICITY> const &);

}}}} // namespace lsst::meas::multifit::grid
