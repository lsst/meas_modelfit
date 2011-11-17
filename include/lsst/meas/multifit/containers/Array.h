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

#ifndef LSST_MEAS_MULTIFIT_CONTAINERS_Array
#define LSST_MEAS_MULTIFIT_CONTAINERS_Array

#include "lsst/meas/multifit/constants.h"
#include "lsst/pex/exceptions.h"
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/format.hpp>
#include <boost/mpl/bool.hpp>
#include <vector>

namespace lsst { namespace meas { namespace multifit { namespace containers {

enum ArrayIndexEnum { NO_INDEX, SORTED, UNSORTED };

namespace detail {

template <typename T>
class ArrayIterator : public boost::iterator_adaptor< 
    ArrayIterator<T>, 
    typename std::vector< boost::shared_ptr<T> >::const_iterator, 
    T, boost::use_default, T &
    > {
    typedef typename ArrayIterator<T>::iterator_adaptor_ Adaptor;
public:

    ArrayIterator() : Adaptor() {}

    ArrayIterator(ArrayIterator const & other) : Adaptor(other.base()) {}

    ArrayIterator(typename std::vector< boost::shared_ptr<T> >::const_iterator const & base) :
        Adaptor(base) {}

    operator boost::shared_ptr<T>() const { return *this->base(); }

private:

    friend class boost::iterator_core_access;

    T & dereference() const { return **this->base(); }

};

template <typename T, typename Derived>
class ArrayBase {
protected:
    typedef std::vector< boost::shared_ptr<T> > Internal;
    Derived const & self() const { return static_cast<Derived const &>(*this); }
public:

    typedef T Value;
    typedef boost::shared_ptr<T> Pointer;
    typedef ArrayIterator<T> Iterator;

    typedef typename boost::remove_const<T>::type value_type;
    typedef T & reference;
    typedef typename boost::add_const<T>::type & const_reference;
    typedef Iterator iterator;
    typedef Iterator const_iterator;
    typedef Pointer pointer;
    typedef typename Internal::difference_type difference_type;
    typedef typename Internal::size_type size_type;

    iterator begin() const { return self().begin(); }
    iterator end() const { return self().end(); }

    reference front() const { return *self().begin(); }
    reference back() const { return *(self().end()-1); }

    size_type size() const { return self().end() - self().begin(); }
    
    bool empty() const { return self().end() == self().begin(); }

};

template <typename Iterator>
Iterator findImpl(ID const id, Iterator iter1, Iterator const end, boost::mpl::true_*) {
    Iterator iter2;
    int count, step;
    count = end - iter1;
    while (count > 0) { // binary search
        iter2 = iter1;
        step = count / 2;
        iter2 += step;
        if (iter2->id < id) {
            iter1 = ++iter2;
            count -= step + 1;
        } else {
            count = step;
        }
    }
    if (iter1->id != id) return end;
    return iter1;
}

template <typename Iterator>
Iterator findImpl(ID const id, Iterator iter, Iterator const end, boost::mpl::false_*) {
    while (iter != end) {
        if (iter->id == id) return iter;
        ++iter;
    }
    return end;
}

template <typename T, typename Derived, bool isSorted>
class IndexedArrayBase : public ArrayBase<T,Derived> {
protected:
    typedef ArrayBase<T,Derived> Super;
public:

    typedef typename Super::value_type value_type;
    typedef typename Super::reference reference;
    typedef typename Super::const_reference const_reference;
    typedef typename Super::iterator iterator;
    typedef typename Super::const_iterator const_iterator;
    typedef typename Super::pointer pointer;
    typedef typename Super::difference_type difference_type;
    typedef typename Super::size_type size_type;
    
    /// @brief Return an iterator to the item with the given ID, or end() if it does not exist.
    iterator find(ID const id) const {
        return findImpl(id, this->begin(), this->end(), (typename boost::mpl::bool_<isSorted>::type *)0);
    }

    /// @brief Return a pointer to the item with the given ID, or an empty pointer if it does not exist.
    pointer get(ID const id) const {
        iterator i = this->find(id);
        if (i == this->end()) return pointer();
        return i;
    }

    /// @brief Return a reference to the item with the given ID, or throw if it does not exist.
    reference operator[](ID const id) const {
        iterator i = this->find(id);
        if (i == this->end()) throw LSST_EXCEPT(
            lsst::pex::exceptions::NotFoundException, 
            (boost::format("Item with ID %d not found.") % id).str()
        );
        return *i;
    }

};

template <typename T, typename Derived, ArrayIndexEnum index> struct SelectBase;

template <typename T, typename Derived>
struct SelectBase<T,Derived,NO_INDEX> {
    typedef ArrayBase<T,Derived> Type;
};

template <typename T, typename Derived>
struct SelectBase<T,Derived,SORTED> {
    typedef IndexedArrayBase<T,Derived,true> Type;
};

template <typename T, typename Derived>
struct SelectBase<T,Derived,UNSORTED> {
    typedef IndexedArrayBase<T,Derived,false> Type;
};

} // namespace detail

template <typename T, ArrayIndexEnum index>
class ArrayView : public detail::SelectBase<T,ArrayView<T,index>,index>::Type {
protected:
    typedef typename detail::SelectBase<T,ArrayView<T,index>,index>::Type Super;
public:

    typedef typename Super::iterator iterator;

    ArrayView() : _first(), _last() {}
    ArrayView(iterator first, iterator last) : _first(first), _last(last) {}

    iterator begin() const { return _first; }
    iterator end() const { return _last; }

protected:

    friend class grid::Initializer;

    void operator=(ArrayView const & other) { _first = other._first; _last = other._last;}

    iterator _first;
    iterator _last;
};

template <typename T, ArrayIndexEnum index>
class Array : public detail::SelectBase<T,Array<T,index>,index>::Type {
    typedef typename detail::SelectBase<T,Array<T,index>,index>::Type Super;
    typedef typename Super::Internal Internal;
public:

    typedef typename Super::iterator iterator;

    Array() : _internal() {}
    explicit Array(Internal & internal) : _internal() { internal.swap(_internal); }

    iterator begin() const { return iterator(_internal.begin()); }
    iterator end() const { return iterator(_internal.end()); }

protected:

    friend class grid::Initializer;

    void operator=(Array const & other) { _internal = other._internal; }

    Internal _internal;
};

}}}} // namespace lsst::meas::multifit::containers

#endif // !LSST_MEAS_MULTIFIT_CONTAINERS_Array
