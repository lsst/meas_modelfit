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

#ifndef LSST_MEAS_MULTIFIT_CONTAINERS_Set
#define LSST_MEAS_MULTIFIT_CONTAINERS_Set

#include "lsst/meas/multifit/constants.h"
#include <boost/format.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <map>
#include <stdexcept>
#include <boost/serialization/map.hpp>

namespace boost {
namespace serialization {
    class access;
}}

namespace lsst { namespace meas { namespace multifit { namespace containers {

namespace detail {

template <typename T>
class SetBase {
protected:
    typedef std::map< ID, boost::shared_ptr<T> > Internal;
    typedef std::pair< ID, boost::shared_ptr<T> > Pair;
public:

    typedef ID Key;
    typedef T Value;
    typedef boost::shared_ptr<T> Pointer;

    class Iterator : public boost::iterator_adaptor< 
        Iterator, typename Internal::const_iterator, Value, boost::use_default, Value &
        > {
        typedef typename Iterator::iterator_adaptor_ Adaptor;
    public:

        Iterator() : Adaptor() {}

        Iterator(Iterator const & other) : Adaptor(other.base()) {}

        Iterator(typename Internal::const_iterator const & base) : Adaptor(base) {}

        operator Pointer() const { return this->base()->second; }

    private:

        friend class boost::iterator_core_access;

        Value & dereference() const { 
            return *this->base()->second;
        }

    };

    typedef Key key_type;
    typedef typename boost::remove_const<Value>::type value_type;
    typedef Value & reference;
    typedef typename boost::add_const<Value>::type & const_reference;
    typedef Iterator iterator;
    typedef iterator const_iterator;
    typedef Pointer pointer;
    typedef typename Internal::difference_type difference_type;
    typedef typename Internal::size_type size_type;

    SetBase() : _internal() {}

    template <typename U>
    SetBase(SetBase<U> const & other) : _internal(other._internal) {}

    iterator begin() const { return iterator(_internal.begin()); }
    iterator end() const { return iterator(_internal.end()); }

    size_type size() const { return _internal.size(); }
    
    bool empty() const { return _internal.empty(); }

    /// @brief Return an iterator to the item with the given ID, or end() if it does not exist.
    iterator find(key_type k) const { return iterator(_internal.find(k)); }

    /// @brief Return a pointer to the item with the given ID, or an empty pointer if it does not exist.
    pointer get(key_type k) const {
        iterator i = find(k);
        if (i == end()) return pointer();
        return i;
    }

    /// @brief Return a reference to the item with the given ID, or throw if it does not exist.
    reference operator[](key_type k) const { 
        iterator i = find(k);
        if (i == end()) throw LSST_EXCEPT(
            lsst::pex::exceptions::NotFoundException, 
            (boost::format("Item with ID %d not found.") % k).str()
        );
        return *i;
    }

protected:

    void operator=(SetBase const & other) {
        _internal = other._internal;
    }

    Internal _internal;
};

} // namespace detail

// Override this for classes with real clone() member functions that will be contained in a Set.
template <typename T>
inline boost::shared_ptr<T> sharedClone(T const & v) {
    return boost::shared_ptr<T>(new T(v));
}

template <typename T>
class ImmutableSet : public detail::SetBase<T> {
protected:
    typedef detail::SetBase<T> Super;
    typedef typename Super::Internal Internal;
    typedef typename Super::Pair Pair;
public:

    typedef typename Super::key_type key_type;
    typedef typename Super::value_type value_type;
    typedef typename Super::reference reference;
    typedef typename Super::const_reference const_reference;
    typedef typename Super::iterator iterator;
    typedef typename Super::const_iterator const_iterator;
    typedef typename Super::Pointer pointer;
    typedef typename Super::difference_type difference_type;
    typedef typename Super::size_type size_type;

    ImmutableSet() : Super() {}

    ImmutableSet(ImmutableSet const & other) : Super(other) {}

    ImmutableSet(Super const & other) : Super(other) {}

};

template <typename T>
class MutableSet : public detail::SetBase<T> {
protected:
    typedef detail::SetBase<T> Super;
    typedef typename Super::Internal Internal;
    typedef typename Super::Pair Pair;

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, unsigned int const version) {
        ar & Super::_internal;
    }

public:

    typedef typename Super::key_type key_type;
    typedef typename Super::value_type value_type;
    typedef typename Super::reference reference;
    typedef typename Super::const_reference const_reference;
    typedef typename Super::iterator iterator;
    typedef typename Super::const_iterator const_iterator;
    typedef typename Super::Pointer pointer;
    typedef typename Super::difference_type difference_type;
    typedef typename Super::size_type size_type;

    MutableSet() : Super() {}

    MutableSet(MutableSet const & other) : Super(other) {}

    MutableSet(Super const & other) : Super(other) {}

    MutableSet & operator=(MutableSet const & other) {
        Super::operator=(other);
        return *this;
    }

    MutableSet & operator=(Super const & other) {
        Super::operator=(other);
        return *this;
    }

    size_type erase(key_type k) { return this->_internal.erase(k); }
    void erase(iterator i) { this->_internal.erase(i.base()); }
    void erase(iterator i1, iterator i2) { 
        this->_internal.erase(i1.base(), i2.base()); 
    }
    void clear() { this->_internal.clear(); }

    std::pair<iterator,bool> insert(const_reference v) {
        return insert(sharedClone(v));
    }

    std::pair<iterator,bool> insert(pointer v) {
        std::pair<typename Internal::iterator,bool> r = this->_internal.insert(Pair(v->id, v));
        return std::pair<iterator,bool>(iterator(r.first), r.second);
    }

    iterator insert(iterator i, const_reference v) { 
        return iterator(this->_internal.insert(i.base(), Pair(v.id, sharedClone(v)))); 
    }

    iterator insert(iterator i, pointer v) { 
        return iterator(this->_internal.insert(i.base(), Pair(v->id, v))); 
    }

    template <typename Iterator_>
    void insert(Iterator_ iter, Iterator_ last) {
        while (iter != last) {
            insert(*iter);
            ++iter;
        }
    }

};

}}}} // namespace lsst::meas::multifit::containers

#endif // !LSST_MEAS_MULTIFIT_CONTAINERS_Set
