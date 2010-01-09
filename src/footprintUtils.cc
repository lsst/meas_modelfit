#include "lsst/meas/multifit/footprintUtils.h"

namespace multifit = lsst::meas::multifit;

lsst::afw::detection::Footprint::Ptr multifit::makeFootprint(
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) {
    lsst::afw::detection::Footprint::Ptr fp(new lsst::afw::detection::Footprint());
    lsst::afw::geom::ellipses::Core::RadialFraction rf(ellipse.getCore());

    lsst::afw::geom::BoxD envelope = ellipse.computeEnvelope();
    int const yEnd = static_cast<int>(envelope.getMaxY()) + 1;
    int const xEnd = static_cast<int>(envelope.getMaxX()) + 1;
    lsst::afw::geom::ExtentD dp(lsst::afw::geom::PointD(0) -ellipse.getCenter());
    for (int y = static_cast<int>(envelope.getMinY()); y<yEnd; ++y) {
        int x = static_cast<int>(envelope.getMinX());
        while (rf(lsst::afw::geom::PointD::make(x,y) + dp) > 1.0) {
            if (x >= xEnd) {
                if (++y >= yEnd) 
                    return fp;
                x = static_cast<int>(envelope.getMinX());
            } else {
                ++x;
            }
        }
        int start = x;
        while (rf(lsst::afw::geom::PointD::make(x,y) + dp) <= 1.0 && x < xEnd) 
            ++x;

        fp->addSpan(y, start, x-1);
    }
    return fp;
}


//explicit template instantiations

template lsst::afw::detection::Footprint::Ptr multifit::clipAndMaskFootprint<lsst::afw::image::MaskPixel>(
    FootprintConstPtr const &,
    lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::Ptr const &
);

template class multifit::CompressFunctor<float>;
template class multifit::CompressFunctor<double>;

template class multifit::ExpandFunctor<float>;
template class multifit::ExpandFunctor<double>;


