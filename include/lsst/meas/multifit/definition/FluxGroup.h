// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_FluxGroup
#define LSST_MEAS_MULTIFIT_DEFINITION_FluxGroup

#include "lsst/meas/multifit/constants.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/serialization/nvp.hpp>

namespace boost {
namespace serialization {
    class access;
}}

namespace lsst { namespace meas { namespace multifit {

namespace grid {
class Initializer;

}

namespace detail {

class FluxGroupBase {
public:

    ID const id;

    /**
     *  @brief Return whether the object is time-variable (has different coefficients
     *         in each exposure).
     */
    bool const isVariable() const { return _variable; };

    /**
     *  @brief Return the maximum value of the norm of the morphology coefficients over the norm of the
     *         flux coefficients.
     */
    double const getMaxMorphologyRatio() const { return _maxMorphologyRatio; }

protected:
    friend class grid::Initializer;
    FluxGroupBase(ID const id_): id(id_) {}
    FluxGroupBase(ID const id_, double maxMorphologyRatio, bool variable) :
        id(id_), _variable(variable), _maxMorphologyRatio(maxMorphologyRatio)
    {}

    FluxGroupBase(FluxGroupBase const & other) :
        id(other.id), _variable(other._variable), _maxMorphologyRatio(other._maxMorphologyRatio)
    {}

    bool _variable;
    double _maxMorphologyRatio;
};

} // namespace detail

namespace definition {

class FluxGroup : public detail::FluxGroupBase {
public:
    
    typedef boost::shared_ptr< FluxGroup > Ptr;
    typedef boost::shared_ptr< FluxGroup const > ConstPtr;

#ifndef SWIG // these are wrapped explicitly; SWIG is confused by the typedefs and "bool &"

    // Use const accessors from base class.
    using detail::FluxGroupBase::isVariable;
    using detail::FluxGroupBase::getMaxMorphologyRatio;

    /**
     *  @brief Return and/or set whether the object is time-variable (has different coefficients
     *         in each exposure).
     */
    bool & isVariable() { return _variable; };

    /**
     *  @brief Return and/or set the maximum value of the norm of the morphology coefficients
     *         over the norm of the flux coefficients.
     */
    double & getMaxMorphologyRatio() { return _maxMorphologyRatio; }


    FluxGroup(FluxGroupBase const &other) : detail::FluxGroupBase(other) {}    
#endif 

    /**
     *  @brief Deep-copy.
     *
     *  @note Not called clone() because it's not virtual, and it doesn't need to be.
     */
    Ptr copy() const { return Ptr(new FluxGroup(*this)); }

    /**
     *  @brief Create a new SharedElement.
     *
     *  Constructors are private to ensure we only get shared_ptrs to these things.
     */
    static Ptr make(ID id, double maxMorphologyRatio, bool variable) {
        return Ptr(new FluxGroup(id, maxMorphologyRatio, variable));
    }

    //@{
    /**
     *  @brief Setters for variability and max morphology ratio.
     *
     *  These aren't necessary in C++ because the non-const getters return references, but they might
     *  be expected, and they're the only way to do things from Python.
     */
    void setVariable(bool variable) { _variable = variable; }
    void setMaxMorphologyRatio(double maxMorphologyRatio) { _maxMorphologyRatio = maxMorphologyRatio; }
    //@}


    FluxGroup(ID id) : detail::FluxGroupBase(id) {}

    FluxGroup(ID id, double maxMorphologyRatio, bool variable) :
        detail::FluxGroupBase(id, maxMorphologyRatio, variable)
    {}

private:
    friend class grid::Initializer;


    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, unsigned int const version) {
        ar & boost::serialization::make_nvp("isVariable", _variable);
        ar & boost::serialization::make_nvp("maxMorphologyRation", _maxMorphologyRatio);
    }
};

}}}} // namespace lsst::meas::multifit::definition

namespace boost { namespace serialization {    
template <class Archive>
inline void save_construct_data(
    Archive & ar, 
    const lsst::meas::multifit::definition::FluxGroup * flux, 
    unsigned int const version
) {
    ar << flux->id;
}
template <class Archive>
inline void load_construct_data(    
    Archive & ar, 
    lsst::meas::multifit::definition::FluxGroup * flux, 
    unsigned int const version
) {
    lsst::meas::multifit::ID id;
    ar >> id;
    ::new(flux) lsst::meas::multifit::definition::FluxGroup(id);
}

}}
#endif // !LSST_MEAS_MULTIFIT_DEFINITION_FluxGroup
