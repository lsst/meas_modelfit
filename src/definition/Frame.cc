#include "lsst/meas/multifit/definition/Frame.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/ndarray/eigen.h"
#include <limits>
#include "Eigen/Core"

namespace lsst { namespace meas { namespace multifit {

namespace detail {

FrameBase::FrameBase(FrameBase const & other, bool deep) :
    id(other.id), _filterId(other._filterId),
    _wcs(other._wcs), _psf(other._psf), _footprint(other._footprint),
    _data(other._data), _weights(other._weights)
{
    if (deep) {
        if (_wcs) _wcs = _wcs->clone();
        if (_psf) _psf = _psf->clone();
        if (_footprint) _footprint = boost::make_shared<Footprint>(*_footprint);
        if (!_data.getData()) _data = ndarray::copy(_data);
        if (!_weights.getData()) _weights = ndarray::copy(_weights);
    }
}

} // namespace detail

namespace definition {

template <typename PixelT>
Frame Frame::make(
    ID const id,
    lsst::afw::image::Exposure<PixelT> const & exposure,
    Footprint::Ptr const & fp,
    lsst::afw::image::MaskPixel const bitmask
) {
    typedef Pixel Pixel;
    typedef lsst::afw::image::MaskedImage<PixelT> MaskedImage;
    
    Footprint::Ptr maskedFp(new Footprint(*fp));
    maskedFp->intersectMask(*exposure.getMaskedImage().getMask(), bitmask);

    //grab portion of exposure's MaskedImage that matches fp
    lsst::afw::geom::BoxI bbox = maskedFp->getBBox();
    MaskedImage mi(
        exposure.getMaskedImage(), 
        bbox,
        lsst::afw::image::PARENT
    );

    //construct flatten version
    lsst::ndarray::Array<Pixel, 1, 1> data = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(maskedFp->getArea())
    );
    lsst::ndarray::Array<Pixel, 1, 1> variance = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(maskedFp->getArea())
    );
    lsst::afw::detection::flattenArray(
        *maskedFp, mi.getImage()->getArray(), data, mi.getXY0()
    );
    lsst::afw::detection::flattenArray(
        *maskedFp, mi.getVariance()->getArray(), variance, mi.getXY0()
    );
    lsst::ndarray::EigenView<Pixel,1,1,Eigen::ArrayXpr> weights(
        lsst::ndarray::allocate(
            lsst::ndarray::makeVector(maskedFp->getArea())
        )
    );
    weights = variance.asEigen<Eigen::ArrayXpr>().sqrt().inverse();

    Frame frame(
        id, 
        exposure.getFilter().getId(), 
        Wcs::Ptr(), 
        exposure.getPsf()->clone(),
        maskedFp,
        data, 
        weights.shallow()
    );

    return frame;
}

template Frame Frame::make<float>(
    ID const, 
    lsst::afw::image::Exposure<float> const &, Footprint::Ptr const &,
    lsst::afw::image::MaskPixel const
);
template Frame Frame::make<double>(
    ID const, 
    lsst::afw::image::Exposure<double> const &, Footprint::Ptr const &, 
    lsst::afw::image::MaskPixel const
);

std::ostream & operator<<(std::ostream & os, Frame const & frame) {
    std::string filterName("undefined");
    try {
        filterName = lsst::afw::image::Filter(frame.getFilterId()).getName();
    } catch (lsst::pex::exceptions::NotFoundException &) {}
    return os << "Frame " << frame.id << " (@" << (&frame) << ") = {" << filterName 
              << ", " << frame.getFootprint()->getArea() << "pix}\n";
}

}}}} // namespace lsst::meas::multifit::definition
