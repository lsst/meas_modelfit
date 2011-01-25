// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_Set
#define LSST_MEAS_MULTIFIT_DEFINITION_Set

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/format.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <stdexcept>

namespace lsst { namespace meas { namespace multifit { namespace definition {

/**
 *  A set/map hybrid that sorts based a member of
 *  the container value_type.
 *
 *  This is mostly just a thin wrapper around a Boost.MultiIndex
 *  container that const_casts references to be mutable (since the
 *  key part is still const).
 */
template <typename Value_, typename Key_, Key_ const Value_::*Member_>
class Set {

    typedef Key_ Key;
    typedef Value_ Value;

    typedef boost::multi_index_container<
        Value,
        boost::multi_index::indexed_by<
            boost::multi_index::ordered_unique<
                boost::multi_index::member< Value, Key const, Member_>
                >
            >
        > Internal;

public:

    typedef typename Internal::value_type value_type;
    typedef Value_ & reference;
    typedef Value_ const & const_reference;
    typedef typename Internal::iterator const_iterator;

    class iterator : 
        public boost::iterator_adaptor< iterator, const_iterator, value_type, boost::use_default, reference>
    {
    public:
        iterator() : iterator::iterator_adaptor_() {}
        iterator(iterator const & other) : iterator::iterator_adaptor_(other) {}

        operator const_iterator() const { return this->base(); }
        
    private:
        friend class boost::iterator_core_access;
        template <typename OtherValue_, typename OtherKey_, OtherKey_ const OtherValue_::*> friend class Set;

        reference dereference() const { return const_cast<reference>(*this->base()); }

        explicit iterator(const_iterator const & other) : iterator::iterator_adaptor_(other) {}
    };

    typedef typename Internal::pointer pointer;
    typedef typename Internal::difference_type difference_type;
    typedef typename Internal::size_type size_type;

    typedef Key key_type;

    Set() : _internal() {}

    template <typename Iterator_>
    Set(Iterator_ first, Iterator_ last) : _internal(first, last) {}

    iterator begin() { return iterator(_internal.begin()); }
    iterator end() { return iterator(_internal.end()); }

    const_iterator begin() const { return _internal.begin(); }
    const_iterator end() const { return _internal.end(); }

    size_type size() const { return _internal.size(); }
    size_type max_size() const { return _internal.max_size(); }
    
    bool empty() const { return _internal.empty(); }

    size_type erase(key_type k) { return _internal.erase(k); }
    void erase(iterator i) { _internal.erase(i.base()); }
    void erase(iterator i1, iterator i2) { _internal.erase(i1.base(), i2.base()); }
    void clear() { _internal.clear(); }

    iterator find(key_type k) { return iterator(_internal.find(k)); }
    const_iterator find(key_type k) const { return _internal.find(k); }

    size_type count(key_type k) const { return _internal.count(k); }

    std::pair<iterator,bool> insert(const_reference v) {
        std::pair<const_iterator,bool> r = _internal.insert(v);
        return std::pair<iterator,bool>(iterator(r.first), r.second);
    }

    template <typename Iterator_>
    void insert(Iterator_ i1, Iterator_ i2) { _internal.insert(i1, i2); }

    iterator insert(iterator i, const_reference v) { return iterator(_internal.insert(i.base(), v)); }

    const_reference operator[](key_type k) const { return *_get(k); }
    reference operator[](key_type k) { return const_cast<reference>(*_get(k)); }

private:

    const_iterator _get(key_type k) const {
        const_iterator i = _internal.find(k);
        if (i == _internal.end()) throw std::invalid_argument(
            boost::str(boost::format("Key %d not found in Set.") % k)
        );
        return i;
    }

    Internal _internal;
};

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Set
