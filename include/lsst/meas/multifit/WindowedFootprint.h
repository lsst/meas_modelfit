#ifndef LSST_MEAS_MULTIFIT_WINDOWED_FOOTPRINT_H
#define LSST_MEAS_MULTIFIT_WINDOWED_FOOTPRINT_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <ndarray_fwd.hpp>

#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/geom/Box.h"

namespace lsst {
namespace meas {
namespace multifit {

class WindowedFootprint {
public:
    typedef boost::shared_ptr<WindowedFootprint> Ptr;
    typedef boost::shared_ptr<const WindowedFootprint> ConstPtr;

    WindowedFootprint(
        lsst::afw::detection::Footprint const & fp, 
        lsst::afw::geom::BoxI const & window
    );

    int const getNpix() const {return _nPix;}

    template <typename T, typename U, int C, int D>
    void compress (
        ndarray::Array<T, 2, C> const & src, 
        ndarray::Array<U, 1, D> const & dest
    ) const {
        if(_window.getWidth() != src.template getSize<1>() ||
            _window.getHeight() != src.template getSize<0>()){
            LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "src must have same dimensions as window"
            );
        }
        if(_nPix != dest.size()) {
            LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "length of dest vector must be equal to number of pixels in footprint"
            );
        }

        Iterator spanIter(_spanList.begin());
        Iterator const & spanEnd(_spanList.end());

        typename ndarray::Array<U,1,C>::Iterator destIter(dest.begin());
        for (; spanIter != spanEnd; ++spanIter) {
            WindowedSpan const & span(**spanIter);
            int spanWidth = span.getWidth();
            if (span.inWindow()) {
                ndarray::Array<T,1,C> row = src[span.getY()];
                std::copy(
                    row.begin() + span.getX0(), 
                    row.begin() + spanWidth,
                    destIter
                );            
            }
            destIter += spanWidth;
        }
    }

    template <typename T, typename U, int N, int C, int D> 
    void compress (
        ndarray::Array<T, N, C> const & src,
        ndarray::Array<U, N - 1, D> const & dest
    ) const {
        if (src.size() != dest.size()) {
            LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Outer dimmension of src and dest do not match");
        }

        typename ndarray::Array<T,N,C>::Iterator const & srcEnd(src.end());
        typename ndarray::Array<T,N,C>::Iterator srcIter(src.begin());
        typename ndarray::Array<U,N-1,C>::Iterator destIter(dest.begin());
        
        for (; srcIter != srcEnd; ++srcIter, ++destIter) {
            compress(*srcIter,*destIter);
        }
    }


    template <typename T, typename U, int C, int D> 
    void expand (
        ndarray::Array<T, 1, C> const & src,
        ndarray::Array<U, 2, D> const & dest
    ) const {
        if(_window.getWidth() != dest.template getSize<1>() ||
            _window.getHeight() != dest.template getSize<0>()){
            LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "dest must have same dimensions as window"
            );
        }
        if(_nPix != src.size()) {
            LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "length of src vector must be equal to number of pixels in footprint"
            );
        }   

        Iterator spanIter(_spanList.begin());
        Iterator const & spanEnd(_spanList.end());

        typename ndarray::Array<T,1,C>::Iterator srcIter(dest.begin());
        for (; spanIter != spanEnd; ++spanIter) {
            WindowedSpan const & span(**spanIter);
            int spanWidth = span.getWidth();
            if (span.inWindow()) {
                ndarray::Array<U,1,C> row = dest[span.getY()];
                std::copy(
                    srcIter, 
                    srcIter + spanWidth,
                    row.begin() + span.getX0()                
                );            
            }
            srcIter += spanWidth;
        }
    }

    template <typename T, typename U, int N, int C, int D> 
    void expand(
        ndarray::Array<T, N - 1, C> const & src,
        ndarray::Array<U, N, D> const & dest
    ) const {
        if (src.size() != dest.size()) {
            LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
            "Outer dimmension of src and dest do not match");
        }
        typename ndarray::Array<T,N-1,C>::Iterator const & srcEnd(src.end());
        typename ndarray::Array<T,N-1,C>::Iterator srcIter(src.begin());
        typename ndarray::Array<U,N,C>::Iterator destIter(dest.begin());

        for (; srcIter != srcEnd; ++srcIter, ++destIter) {
            expand(*srcIter,*destIter);
        }
    }



private:
    typedef lsst::afw::detection::Span Span;

    class WindowedSpan : public Span {
    public:
        typedef boost::shared_ptr<WindowedSpan> Ptr;

        WindowedSpan(int const & y, int const & x0, int const & x1, bool inWindow)
            : Span(y, x0, x1), _inWindow(inWindow) {}
        
        
        bool const & inWindow() const {return _inWindow;}
    private:
        bool _inWindow;
    };

    typedef lsst::afw::detection::Footprint Footprint;

    typedef std::list<WindowedSpan::Ptr> SpanMap;
    typedef SpanMap::const_iterator Iterator;

    int _nPix;
    SpanMap _spanList;
    lsst::afw::geom::BoxI _window;    
};

}}}
#endif
