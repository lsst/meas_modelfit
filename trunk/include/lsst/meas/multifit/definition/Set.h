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
 *  A never-const shared_ptr-based set/map hybrid that sorts based the id member of Value.
 */
template <typename Value_>
class Set {
    typedef std::map< ID, boost::shared_ptr<Value_> > Internal;
    typedef std::pair< ID, boost::shared_ptr<Value_> > Pair;
public:

    typedef ID Key;
    typedef Value_ Value;
    typedef boost::shared_ptr<Value_> Pointer;

    class Iterator : public boost::iterator_adaptor< 
        Iterator, typename Internal::const_iterator, Value, boost::use_default, Value &
    > {
        typedef typename Iterator::iterator_adaptor_ Adaptor;
    public:
        Iterator() : Adaptor() {}

        Iterator(Iterator const & other) : Adaptor(other.base()) {}

        operator Pointer() const { return this->base()->second; }

    private:

        friend class boost::iterator_core_access;

        template <typename OtherValue_> friend class Set;

        Value & dereference() const { 
            return *this->base()->second;
        }

        Iterator(typename Internal::const_iterator const & adapted) : Adaptor(adapted) {}
    };

    typedef Key key_type;
    typedef Value value_type;
    typedef Value & reference;
    typedef Value const & const_reference;
    typedef Iterator iterator;
    typedef Iterator const_iterator;
    typedef Pointer pointer;
    typedef typename Internal::difference_type difference_type;
    typedef typename Internal::size_type size_type;

    Set() : _internal() {}

    template <typename Iterator_>
    Set(Iterator_ iter, Iterator_ last) : _internal() { insert(iter, last); }

    iterator begin() const { return iterator(_internal.begin()); }
    iterator end() const { return iterator(_internal.end()); }

    size_type size() const { return _internal.size(); }
    size_type max_size() const { return _internal.max_size(); }
    
    bool empty() const { return _internal.empty(); }

    size_type erase(key_type k) { return _internal.erase(k); }
    void erase(iterator i) { _internal.erase(i.base()); }
    void erase(iterator i1, iterator i2) { 
        _internal.erase(i1.base(), i2.base()); 
    }
    void clear() { _internal.clear(); }

    iterator find(key_type k) const { return iterator(_internal.find(k)); }

    size_type count(key_type k) const { return _internal.count(k); }

    std::pair<iterator,bool> insert(const_reference v) {
        std::pair<typename Internal::iterator,bool> r = _internal.insert(Pair(v.id, Pointer(new Value(v))));
        return std::pair<iterator,bool>(iterator(r.first), r.second);
    }

    std::pair<iterator,bool> insert(pointer v) {
        std::pair<typename Internal::iterator,bool> r = _internal.insert(Pair(v->id, v));
        return std::pair<iterator,bool>(iterator(r.first), r.second);
    }

    iterator insert(iterator i, const_reference v) { 
        return iterator(_internal.insert(i.base(), Pair(v.id, Pointer(new Value(v))))); 
    }

    iterator insert(iterator i, pointer v) { 
        return iterator(_internal.insert(i.base(), Pair(v->id, v))); 
    }

    template <typename Iterator_>
    void insert(Iterator_ iter, Iterator_ last) {
        while (iter != last) {
            insert(*iter);
            ++iter;
        }
    }

    reference operator[](key_type k) const { return *_get(k); }

private:

    iterator _get(key_type k) const {
        typename Internal::const_iterator i = _internal.find(k);
        if (i == _internal.end()) throw std::invalid_argument(
            boost::str(boost::format("Key %d not found in Set.") % k)
        );
        return iterator(i);
    }

    Internal _internal;
};

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Set
