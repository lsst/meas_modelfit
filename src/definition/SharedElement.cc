#include "lsst/meas/multifit/definition/SharedElement.h"

namespace lsst { namespace meas { namespace multifit { namespace definition {

template <SharedElementType E>
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & component) {
    static char const * names[] = { "Position", "Radius", "Ellipticity" };
    os << names[E] << " (@" << (&component) << ") = ";
    detail::SharedElementTraits<E>::printValue(os, component.getValue());
    return os << (component.isActive() ? " (active) " : " (inactive) ") 
              << "(" << component.getBounds() << ")";
}

template std::ostream & operator<<(std::ostream &, SharedElement<POSITION> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<RADIUS> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<ELLIPTICITY> const &);

}}}} // namespace lsst::meas::multifit::definition
