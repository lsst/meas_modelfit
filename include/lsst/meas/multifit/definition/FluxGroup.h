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

namespace lsst { namespace meas { namespace multifit {

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
    
    FluxGroupBase(ID const id_, double maxMorphologyRatio, bool isVariable) :
        id(id_), _variable(isVariable), _maxMorphologyRatio(maxMorphologyRatio)
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
    static Ptr make(ID id, double maxMorphologyRatio, bool isVariable) {
        return Ptr(new FluxGroup(id, maxMorphologyRatio, isVariable));
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

private:

    FluxGroup(FluxGroup const & other) : detail::FluxGroupBase(other) {}
    
    FluxGroup(ID id, double maxMorphologyRatio, bool isVariable) :
        detail::FluxGroupBase(id, isVariable, maxMorphologyRatio)
    {}

};

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_FluxGroup
