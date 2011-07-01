#include "lsst/meas/multifit/definition/Definition.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/ndarray/eigen.h"
#include <limits>
#include <Eigen/Array>

namespace lsst { namespace meas { namespace multifit { namespace definition {

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
