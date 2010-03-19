// -*- lsst-c++ -*-
#include <lsst/meas/multifit/WindowedFootprint.h>
#include <ndarray.hpp>
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;
namespace afwDet = lsst::afw::detection;
/**
 * Construct a WindowedFootprint
 *
 * The construction results in a deep copy of the span's in the
 * Footprint, and as necessary, splits any segments which intersect the bounding
 * box. Every span is then labeled as being inside or outside the bounding box.
 * 
 * @note The resulting WindowedFootprint's pixel indexes are relative to the 
 *      window. So given a window spanning from (4,4) to (10,10) inclusive, and
 *      fp containing the pixel indexes (4,4) and (10, 10), the resulting 
 *      WindowedFootprint, will contain the pixel index (0,0) and (6,6)
 *
 * @param fp The original, footprint to constrain
 * @param window The box defining the bounds of the desired work area
 */
multifit::WindowedFootprint::WindowedFootprint (
    lsst::afw::detection::Footprint const & fp,
    lsst::afw::geom::BoxI const & window
) : 
    _nWindowedPix(0),
    _nPix(fp.getNpix()), 
    _window(lsst::afw::geom::PointI(0), window.getDimensions()),
    _spanList(new SpanList())
{
    afwDet::Footprint::SpanList spanList = fp.getSpans();
    // Define the minimum relative to the input window    
    lsst::afw::geom::Point2I min = window.getMin();
    // Define the maximum relative to 0     
    lsst::afw::geom::Point2I max = _window.getMax();

    for (afwDet::Footprint::SpanList::const_iterator i(spanList.begin()), 
         end(spanList.end()); i != end; ++i
    ) {
        int x0 = (*i)->getX0() - min.getX();
        int x1 = (*i)->getX1() - min.getX();
        int y = (*i)->getY() - min.getY();
        if (y < 0 || y > max.getY() || x1 < 0 || x0 > max.getX()) {
            //The span is completely outside the window
            _spanList->push_back(WindowedSpan(y, x0, x1, false));
            
            //nothing more to be done for this span
            continue;
        }
        
        if (x0 < 0) {
            // Span extends to the left of the window.
            // insert just the portion of the span extending to the left
            _spanList->push_back(WindowedSpan(y, x0, -1, false));
           
            //push x0 to start of window
            x0 = 0;
        }
        if (x1 > max.getX()) {
            //span extends to the right of the window
            
            //insert portion within the window            
            _spanList->push_back(WindowedSpan(y, x0, max.getX(), true));
            _nWindowedPix += max.getX() - x0 +1;
    
            //shift start to just to the right of the window
            x0 = max.getX() + 1;

            //insert the portion of the span that is outside the window
            _spanList->push_back(WindowedSpan(y, x0, x1, false));
            
            x1 = max.getX();
        }
        else if(x0 <= x1) {
            //insert portion within the window            
            _spanList->push_back(WindowedSpan(y, x0, x1, true));
            _nWindowedPix += x1 - x0 +1;
        }
    }
}
