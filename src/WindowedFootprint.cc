#include <lsst/meas/multifit/WindowedFootprint.h>
#include <ndarray.hpp>
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::WindowedFootprint::WindowedFootprint (
    lsst::afw::detection::Footprint const & fp,
    lsst::afw::geom::BoxI const & window
) : 
    _nPix(0), 
    _window(lsst::afw::geom::PointI(0), window.getDimensions()) 
{
    Footprint::SpanList spanList = fp.getSpans();
    Footprint::SpanList::const_iterator i(spanList.begin());
    Footprint::SpanList::const_iterator const & end(spanList.end());
    
    lsst::afw::geom::Point2I min = window.getMin();
    lsst::afw::geom::Point2I max = _window.getMax();
    WindowedSpan::Ptr toAdd;
    int x0, x1, y;
    for (; i != end; ++i) {
        x0 = (*i)->getX0() - min.getX();
        x1 = (*i)->getX1() - min.getX();
        y = (*i)->getY() - min.getY();
        if (y < 0 || y > max.getY() || x1 < 0 || x0 > max.getX()) {
            //The span is completely outside the window
            toAdd = boost::make_shared<WindowedSpan>(y, x0, x1, false);
            _spanList.push_back(toAdd);
            
            //nothing more to be done for this span
            continue;
        }
        
        if (x0 < 0) {
            // Span extends to the left of the window.
            // insert just the portion of the span extending to the left
            toAdd = boost::make_shared<WindowedSpan>(y, x0, -1, false);
            _spanList.push_back(toAdd);
           
            //push x0 to start of window
            x0 = 0;
        }

        if (x1 > max.getX()) {
            //span extends to the right of the window
            //insert the portion of the span that is within the window
            toAdd = boost::make_shared<WindowedSpan>(y, x0, max.getX(), true);
            _spanList.push_back(toAdd);

            _nPix += toAdd->getWidth();     

            //insert the portion of the span extending to the right
            toAdd = boost::make_shared<WindowedSpan>(y, max.getX() +1, x1, false);
            _spanList.push_back(toAdd);
        } else {
            //insert portion within the window            
            toAdd = boost::make_shared<WindowedSpan>(y, x0, x1, true);
            _spanList.push_back(toAdd);

            _nPix += toAdd->getWidth();
        }
    }
}
