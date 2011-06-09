#include "lsst/meas/multifit/definition/Definition.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/ndarray/eigen.h"
#include <limits>
#include <Eigen/Array>

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
        // TODO: copy footprint when copy ctor becomes available in afw
        if (!_data.getData()) _data = ndarray::copy(_data);
        if (!_weights.getData()) _weights = ndarray::copy(_weights);
    }
}

} // namespace detail


namespace definition {

template <SharedElementType E>
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & element) {
    static char const * names[] = { "Position", "Radius", "Ellipticity" };
    os << names[E] << " (@" << (&element) << ") = ";
    detail::SharedElementTraits<E>::printValue(os, element.getValue());
    return os << (element.isActive() ? " (active) " : " (inactive) ") 
              << "(" << element.getBounds() << ")";
}

template std::ostream & operator<<(std::ostream &, SharedElement<POSITION> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<RADIUS> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<ELLIPTICITY> const &);

ObjectComponent ObjectComponent::makeStar(
    ID const id, 
    lsst::afw::geom::Point2D const & position, 
    bool const isVariable,
    bool const isPositionActive
) {
    ObjectComponent r(id);
    r.getPositionElement() = PositionElement::make(position, isPositionActive);
    r.isVariable() = isVariable;
    return r;
}

ObjectComponent ObjectComponent::makeGalaxy(
    ID const id,
    ModelBasis::Ptr const & basis,
    lsst::afw::geom::ellipses::Ellipse const & ellipse,
    bool const isEllipticityActive,
    bool const isRadiusActive,
    bool const isPositionActive
) {
    ObjectComponent r(id);
    EllipseCore core(ellipse.getCore());
    r.getPositionElement() = PositionElement::make(ellipse.getCenter(), isPositionActive);
    r.getEllipticityElement() = EllipticityElement::make(core.getEllipticity(), isEllipticityActive);
    r.getRadiusElement() = RadiusElement::make(core.getRadius(), isRadiusActive);
    r.getBasis() = basis;
    return r;
}

std::ostream & operator<<(std::ostream & os, ObjectComponent const & obj) {
    os << "ObjectComponent " << obj.id << "(@" << (&obj) << ") = {"
       << (obj.isVariable() ? "variable" : "nonvariable") << ", Rx" << obj.getRadiusFactor() << "}:\n";
    if (obj.getPositionElement()) os << "    " << (*obj.getPositionElement()) << "\n";
    if (obj.getRadiusElement()) 
        os << "    " << (*obj.getRadiusElement()) << " x " << obj.getRadiusFactor() << "\n";
    if (obj.getEllipticityElement()) os << "    " << (*obj.getEllipticityElement()) << "\n";
    return os;
}

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
    lsst::ndarray::EigenView<Pixel, 1, 1> weights(
        lsst::ndarray::allocate(
            lsst::ndarray::makeVector(maskedFp->getArea())
        )
    );
    weights = lsst::ndarray::viewAsEigen(variance).cwise().sqrt().cwise().inverse();
    Frame frame(
        id, 
        exposure.getFilter().getId(), 
        Wcs::Ptr(), 
        exposure.getPsf()->clone(),
        maskedFp,
        data, 
        weights.getArray()
    );

    return frame;
}

template Frame Frame::make<float>(
    ID const, lsst::afw::image::Exposure<float> const &, Footprint::Ptr const &,
    lsst::afw::image::MaskPixel const
);
template Frame Frame::make<double>(
    ID const, lsst::afw::image::Exposure<double> const &, Footprint::Ptr const &, 
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

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), _wcs(other._wcs) {}

}}}} // namespace lsst::meas::multifit::definition
