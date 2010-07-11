// -*- lsst-c++ -*-

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
 
/**
 * @file 
 *
 * Declaration of the WindowedFootprint class
 */
#ifndef LSST_MEAS_MULTIFIT_WINDOWED_FOOTPRINT_H
#define LSST_MEAS_MULTIFIT_WINDOWED_FOOTPRINT_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <ndarray_fwd.hpp>
#include <algorithm>

#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/geom/Box.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 * A footprint with a constrained working area
 *
 * The WindowedFootprint is constructed from a Footprint and a BoxI
 * The construction results in a deep copy of the span's in the
 * Footprint, and as necessary, splits any segments which intersect the bounding
 * box. Every span is then labeled as being inside or outside the bounding box.
 *
 * This caching of which pixels of the footprint are within a desired work
 * area, is a small speed optimization over iterating over a conventional
 * footprint, and having to test at each point whether the current pixel is
 * outside the desired bounding box
 *
 * Within multifit the WindowedFootprint is used to compress or expand arrays
 * according to the layout of a footprint. Given a 2d array (image or matrix),
 * the Windowed footprint creates a 1d array, where only those points covered
 * by the footprint are present. The inverse operation is also supported, 
 * expanding a compressed array to a which include all pixels within the window 
 * of the footprint, and all the pixels covered by the footprint are set to the
 * value of the corresponding pixel in the compressed array.
 */
class WindowedFootprint {
public:
    typedef boost::shared_ptr<WindowedFootprint> Ptr;
    typedef boost::shared_ptr<const WindowedFootprint> ConstPtr;
    
    WindowedFootprint(
        lsst::afw::detection::Footprint const & fp, 
        lsst::afw::geom::BoxI const & window
    );

    /**
     * Number of pixels in the entire footprint
     *
     * The naming convention is analogous to lsst::afw::detection::Footprint
     */
    int const &getNpix() const {return _nPix;}

    /**
     * Number of pixels inside the bounds of the window     
     */
    int const & getNwindowedPix() const {return _nWindowedPix;} 

    /**
     * The specified bounding box
     */
    lsst::afw::geom::BoxI const & getWindow() const {return _window;}

    /**
     * Footprint-Compress a 2-d array into a 1-d array
     *
     * After this function is called, dest will contain all values of src
     * covered by this WindowedFootprint, flattened into a 1-d array
     *
     * @param src must be externally allocated with dimensions (window.getHeight(), 
     *      window.getWidth())
     * @param dest must be externally allocated with size (nPix)
     * @throw lsst::pex::exceptions::LengthErrorException if src or dest do
     *      not have the right dimensions
     */
    template <typename T, typename U, int C>
    void compress (
        ndarray::Array<T, 2, C> const & src, 
        ndarray::Array<U, 1, 1> const & dest
    ) const {
        if(_window.getWidth() != src.template getSize<1>() ||
            _window.getHeight() != src.template getSize<0>()){
            LSST_EXCEPT(
                lsst::pex::exceptions::LengthErrorException,
                "src must have same dimensions as window"
            );
        }
        if(_nPix != dest.size()) {
            LSST_EXCEPT(
                lsst::pex::exceptions::LengthErrorException,
                "length of dest vector must be equal to number of pixels in footprint"
            );
        }

        typename ndarray::Array<U, 1, 1>::Iterator destIter = dest.begin();
        typename ndarray::Array<T, 2, C>::Iterator srcIter = src.begin();
        int lastY = 0; 
        for (SpanIterator i(_spanList->begin()), end(_spanList->end()); i != end; ++i) {
            WindowedSpan const & span(*i);
            int spanWidth = span.getWidth();

            srcIter += span.getY() - lastY;
            lastY = span.getY();

            if (span.inWindow()) {
                typename ndarray::Array<T, 2, C>::Reference::Iterator data(
                    srcIter->begin() + span.getX0()
                );
                std::copy(data, data + spanWidth, destIter);            
            } else {
                std::fill_n(destIter, spanWidth, 0);
            }
            
            destIter += spanWidth;
        }
    }

    /**
     * Footprint-Compress a n-d array into a (n-1)-d array
     *
     * After this function is called, dest will contain all values of src
     * covered by this WindowedFootprint in each 2-d array it contains, 
     * flattened into a (n-1)-d array. 
     *
     * @param src must be externally allocated where the two inner dimensions 
     *      are (window.getHeight(), window.getWidth())
     * @param dest must be externally allocated where the innermost dimension
     *      is nPix, and all outer dimensions are equivalent in size to the 
     *      matching src dimension
     * @throw lsst::pex::exceptions::LengthErrorException if src or dest 
     *      do not have the right dimensions
     */
    template <typename T, typename U, int N, int C, int D> 
    void compress (
        ndarray::Array<T, N, C> const & src,
        ndarray::Array<U, N - 1, D> const & dest
    ) const {
        if (src.size() != dest.size()) {
            LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Outer dimension of src and dest do not match");
        }

        typename ndarray::Array<T,N,C>::Iterator const & srcEnd(src.end());
        typename ndarray::Array<T,N,C>::Iterator srcIter(src.begin());
        typename ndarray::Array<U,N-1,D>::Iterator destIter(dest.begin());
       
        //loop over outer dimension        
        for (; srcIter != srcEnd; ++srcIter, ++destIter) {
            //recurse into inner dimensions
            compress(*srcIter,*destIter);
        }
    }

    /**
     * Footprint-Expand a 1-d array into a 2-d array
     *
     * After this function is called, all values of dest covered by this 
     * WindowedFootprint will be set to the value specified in src, where this
     * footprint specifies the mapping.    
     *
     * @param src must be externally allocated with size (nPix)
     * @param dest must be externally allocated with dimensions 
     *      (window.getHeight(), window.getWidth())
     * @throw lsst::pex::exceptions::LengthErrorException if src or dest 
     *      do not have the right dimensions
     */
    template <typename T, typename U, int C> 
    void expand (
        ndarray::Array<T, 1, 1> const & src,
        ndarray::Array<U, 2, C> const & dest
    ) const {
        if(_window.getWidth() != dest.template getSize<1>() ||
            _window.getHeight() != dest.template getSize<0>()){
            LSST_EXCEPT(
                lsst::pex::exceptions::LengthErrorException,
                "dest must have same dimensions as window"
            );
        }
        if(_nPix != src.size()) {
            LSST_EXCEPT(
                lsst::pex::exceptions::LengthErrorException,
                "length of src vector must be equal to number of pixels in footprint"
            );
        }   

        typename ndarray::Array<T, 1, 1>::Iterator srcIter(src.begin());    
        typename ndarray::Array<U, 2, C>::Iterator destRow(dest.begin());
        int lastY = 0;
        for (SpanIterator i(_spanList->begin()), end(_spanList->end()); i != end; ++i) {
            WindowedSpan const & span(*i);
            int spanWidth = span.getWidth();

            destRow += span.getY() - lastY;
            lastY = span.getY();

            typename ndarray::Array<T, 2, C>::Reference::Iterator data(
                destRow->begin() + span.getX0()
            );

            if (span.inWindow()) {
                std::copy(
                    srcIter, 
                    srcIter + spanWidth,
                    data
                );            
            } else {
                std::fill_n(data, spanWidth, 0);
            }
            srcIter += spanWidth;
        }
    }

    /**
     * Footprint-Expand a (n-1)-d array into a (n)-d array
     *
     * After this function is called, all indexes of each 2-d array in dest 
     * covered by this WindowedFootprint will be set to the value specified in
     * src, where this footprint specifies the mapping. 
     *
     * @param src must be externally allocated where the two inner dimensions 
     *      are (window.getHeight(), window.getWidth())
     * @param dest must be externally allocated where the innermost dimension 
     *      has size nPix, and all outer dimensions are equivalent in size to 
     *      the corresponding src dimension
     * @throw lsst::pex::exceptions::InvalidParameterException if src or dest do
     *     not have the right dimensions
     */
    template <typename T, typename U, int N, int C, int D> 
    void expand(
        ndarray::Array<T, N - 1, C> const & src,
        ndarray::Array<U, N, D> const & dest
    ) const {
        if (src.size() != dest.size()) {
            LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
            "Outer dimension of src and dest do not match");
        }
        typename ndarray::Array<T,N-1,C>::Iterator const & srcEnd(src.end());
        typename ndarray::Array<T,N-1,C>::Iterator srcIter(src.begin());
        typename ndarray::Array<U,N,D>::Iterator destIter(dest.begin());

        //loop over outer dimension
        for (; srcIter != srcEnd; ++srcIter, ++destIter) {
            //recurse to inner dimension
            expand(*srcIter,*destIter);
        }
    }



private:
    typedef lsst::afw::detection::Span Span;

    /**
     * Internal representation of a span which knows whether it lies inside or
     * outside the window of the WindowedFootprint it belongs to
     */
    class WindowedSpan : public Span {
    public:
        typedef boost::shared_ptr<WindowedSpan> Ptr;

        /** 
         * Construct a set of pixel indexes spanning from (x0, y) to (x1, y)
         * inclusive, which state whether the pixels are inside or outside the
         * window of the WindowedFootprint it belongs to
         *
         * @param y row index of the span
         * @param x0 starting column index of the span
         * @param x1 ending column index (inclusive) of the span
         * @param inWindow define whether the span is considered in or out of
         *      bounds
         */
        WindowedSpan(int const & y, int const & x0, int const & x1, bool inWindow)
            : Span(y, x0, x1), _inWindow(inWindow) {}
        
        /**
         * state whether the span is considered in or out of the bounds of the
         * WindowedFootprint it belongs to
         */
        bool const & inWindow() const {return _inWindow;}
    private:
        bool _inWindow;
    };

    typedef lsst::afw::detection::Footprint Footprint;

    typedef std::list<WindowedSpan> SpanList;
    typedef boost::shared_ptr<SpanList> SpanListPtr;
    typedef SpanList::const_iterator SpanIterator;

    int _nWindowedPix, _nPix;
    lsst::afw::geom::BoxI _window;    
    SpanListPtr _spanList;
};

}}}
#endif
