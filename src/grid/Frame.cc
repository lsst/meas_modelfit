#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

Frame::Frame(definition::Frame const & def, int pixelOffset, int filterIndex, int frameIndex) :
    detail::FrameBase(def, true), _pixelOffset(pixelOffset), _pixelCount(_footprint->getNpix()), 
    _filterIndex(filterIndex), _frameIndex(frameIndex)
{}

void Frame::applyWeights(ndarray::Array<Pixel,2,1> const & matrix) const {
    if (!_weights.empty()) {
        matrix.deep() *= _weights;
    }
}

void Frame::applyWeights(ndarray::Array<Pixel,1,0> const & vector) const {
    if (!_weights.empty()) {
        vector.deep() *= _weights;
    }
}

std::ostream & operator<<(std::ostream & os, Frame const & frame) {
    std::string filterName("undefined");
    try {
        filterName = lsst::afw::image::Filter(frame.getFilterId()).getName();
    } catch (lsst::pex::exceptions::NotFoundException &) {}
    return os << "Frame " << frame.id << " (@" << (&frame) << ") = {" << filterName 
              << ", " << frame.getFootprint()->getArea() << "pix}\n";
}

}}}} // namespace lsst::meas::multifit::grid
