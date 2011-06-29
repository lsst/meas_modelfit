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
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & component) {
    static char const * names[] = { "Position", "Radius", "Ellipticity" };
    os << names[E] << " (@" << (&component) << ") = ";
    detail::SharedElementTraits<E>::printValue(os, component.getValue());
    return os << (component.isActive() ? " (active) " : " (inactive) ") 
              << "(" << component.getBounds() << ")";
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
    r.getPosition() = PositionElement::make(position, isPositionActive);
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
    r.getPosition() = PositionElement::make(ellipse.getCenter(), isPositionActive);
    r.getEllipticity() = EllipticityElement::make(core.getEllipticity(), isEllipticityActive);
    r.getRadius() = RadiusElement::make(core.getRadius(), isRadiusActive);
    r.getBasis() = basis;
    return r;
}

std::ostream & operator<<(std::ostream & os, ObjectComponent const & obj) {
    os << "ObjectComponent " << obj.id << "(@" << (&obj) << ") = {"
       << (obj.isVariable() ? "variable" : "nonvariable") << ", Rx" << obj.getRadiusFactor() << "}:\n";
    if (obj.getPosition()) os << "    " << (*obj.getPosition()) << "\n";
    if (obj.getRadius()) os << "    " << (*obj.getRadius()) << " x " << obj.getRadiusFactor() << "\n";
    if (obj.getEllipticity()) os << "    " << (*obj.getEllipticity()) << "\n";
    return os;
}

template <typename PixelT>
Frame Frame::make(
    ID const id,
    lsst::afw::image::Exposure<PixelT> const & exposure,
    Footprint::Ptr const & fp,
    lsst::afw::image::MaskPixel const bitmask,
    bool const usePixelWeights
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
    if(usePixelWeights) {
        weights = lsst::ndarray::viewAsEigen(variance).cwise().sqrt().cwise().inverse();
    } else if (weights.size() > 0) {
        weights.setOnes();
    }

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
    ID const, 
    lsst::afw::image::Exposure<float> const &, Footprint::Ptr const &,
    lsst::afw::image::MaskPixel const,
    bool const
);
template Frame Frame::make<double>(
    ID const, 
    lsst::afw::image::Exposure<double> const &, Footprint::Ptr const &, 
    lsst::afw::image::MaskPixel const,
    bool const
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

template <typename PixelT>
Definition Definition::make(
    afw::image::Exposure<PixelT> const & exposure,
    Footprint::Ptr const & fp,
    afw::geom::Point2D const & position,
    bool const isVariable, 
    bool const isPositionActive,
    lsst::afw::image::MaskPixel const bitmask,
    bool const usePixelWeights
) {    
    //make a point source definition   
    Definition psDefinition;
    psDefinition.frames.insert(
        Frame::make<PixelT>(0, exposure, fp, bitmask, usePixelWeights)
    );
    psDefinition.objects.insert(
        definition::ObjectComponent::makeStar(0, position, isVariable, isPositionActive)
    );
    return psDefinition;
}

template Definition Definition::make<float>(
    lsst::afw::image::Exposure<float> const &,
    lsst::meas::multifit::Footprint::Ptr const &,
    afw::geom::Point2D const&,
    bool const, bool const,
    lsst::afw::image::MaskPixel const,
    bool const
);
template Definition Definition::make<double>(
    lsst::afw::image::Exposure<double> const &,
    lsst::meas::multifit::Footprint::Ptr const &,
    afw::geom::Point2D const&,
    bool const, bool const,
    lsst::afw::image::MaskPixel const,
    bool const
);

template <typename PixelT>
Definition Definition::make(
    afw::image::Exposure<PixelT> const & exposure,
    Footprint::Ptr const & fp,
    ModelBasis::Ptr const & basis,
    afw::geom::ellipses::Ellipse const & ellipse,
    bool const isEllipticityActive,
    bool const isRadiusActive,
    bool const isPositionActive,
    lsst::afw::image::MaskPixel const bitmask,
    bool const usePixelWeights
) {
    //make a single-galaxy definition    
    Definition sgDefinition;
    sgDefinition.frames.insert(
        Frame::make<PixelT>(0, exposure, fp, bitmask, usePixelWeights)
    );
    sgDefinition.objects.insert(
        definition::ObjectComponent::makeGalaxy(
            0, basis, ellipse,
            isEllipticityActive,
            isRadiusActive,
            isPositionActive
        )
    );
    return sgDefinition;
}

template Definition Definition::make<float>(
    lsst::afw::image::Exposure<float> const &,
    lsst::meas::multifit::Footprint::Ptr const &,
    ModelBasis::Ptr const &,
    afw::geom::ellipses::Ellipse const &,
    bool const, bool const, bool const,
    lsst::afw::image::MaskPixel const,
    bool const
);
template Definition Definition::make<double>(
    lsst::afw::image::Exposure<double> const &,
    lsst::meas::multifit::Footprint::Ptr const &,
    ModelBasis::Ptr const &,
    afw::geom::ellipses::Ellipse const &,
    bool const, bool const, bool const,
    lsst::afw::image::MaskPixel const,
    bool const
);

}}}} // namespace lsst::meas::multifit::definition
