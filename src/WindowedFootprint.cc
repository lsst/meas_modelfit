#include <lsst/meas/multifit/WindowedFootprint.h>
#include <ndarray.hpp>
#include <lsst/afw/image/Utils.h>
#include <lsst/pex/exceptions/Runtime.h>

namespace multifit = lsst::meas::multifit;

multifit::WindowedFootprint::WindowedFootprint (
    lsst::afw::detection::Footprint const & fp,
    lsst::afw::image::BBox const & window
) : 
    _nPix(fp.getNpix()), 
    _window(window) 
{
    Footprint::SpanList spanList = fp.getSpans();
    Footprint::SpanList::const_iterator i(spanList.begin());
    Footprint::SpanList::const_iterator const & end(spanList.end());
    
    lsst::afw::image::PointI llc = window.getLLC();
    lsst::afw::image::PointI urc = window.getURC() - llc;   
    int x0, x1, y;
    for (; i != end; ++i) {
        x0 = (*i)->getX0() - llc.getX();
        x1 = (*i)->getX1() - llc.getX();
        y = (*i)->getY() - llc.getY();
        if (y < 0 || y > urc.getY() || x1 < 0 || x0 > urc.getX()) {
            //The span is completely outside the window
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, x1, false)
            ); 

            //nothing more to be done for this span
            continue;
        }
        

        if (x0 < llc.getX()) {
            // Span extends to the left of the window.
            // insert just the portion of the span extending to the left
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, -1, false)
            );
            x0 = 0;
        }

        if (x1 > urc.getX()) {
            //span extends to the right of the window
            //insert the portion of the span that is within the window
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, urc.getX(), true)
            );
                
            //insert the portion of the span extending to the right
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, urc.getX() +1, x1, false)
            );
        } else {
            //insert portion within the window            
            _spanList.push_back(
                boost::make_shared<WindowedSpan>(y, x0, x1, true)
            );
        }
    }
}



