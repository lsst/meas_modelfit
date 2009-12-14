#include <lsst/meas/multifit/WindowedFootprint.h>
#include <ndarray.hpp>
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::WindowedFootprint::WindowedFootprint (
    lsst::afw::detection::Footprint const & fp,
    lsst::afw::geom::Box2I const & window
) : 
    _nPix(fp.getNpix()), 
    _window(window) 
{
    Footprint::SpanList spanList = fp.getSpans();
    Footprint::SpanList::const_iterator i(spanList.begin());
    Footprint::SpanList::const_iterator const & end(spanList.end());
    
    lsst::afw::geom::Point2I min = window.getMin();
    lsst::afw::geom::Extent2I dimension = window.getDimensions();
    int x0, x1, y;
    for (; i != end; ++i) {
        x0 = (*i)->getX0() - min.getX();
        x1 = (*i)->getX1() - min.getX();
        y = (*i)->getY() - min.getY();
        if (y < 0 || y > dimension.getY() || x1 < 0 || x0 > dimension.getX()) {
            //The span is completely outside the window
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, x1, false)
            ); 

            //nothing more to be done for this span
            continue;
        }
        

        if (x0 < min.getX()) {
            // Span extends to the left of the window.
            // insert just the portion of the span extending to the left
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, -1, false)
            );
            x0 = 0;
        }

        if (x1 > dimension.getX()) {
            //span extends to the right of the window
            //insert the portion of the span that is within the window
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, dimension.getX(), true)
            );
                
            //insert the portion of the span extending to the right
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, dimension.getX() +1, x1, false)
            );
        } else {
            //insert portion within the window            
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, x1, true)
            );
        }
    }
}
