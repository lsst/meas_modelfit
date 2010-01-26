// -*- lsst-c++ -*-
/**
 * @file
 * Support for an Exposure characterized by a PSF
 *
 */
#ifndef LSST_MEAS_MULTIFIT_CHARACTERIZED_EXPOSURE_H
#define LSST_MEAS_MULTIFIT_CHARACTERIZED_EXPOSURE_H

#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/PSF.h"

namespace lsst { namespace meas { namespace multifit {

/**
 * An lsst exposure characterized by a PSF
 *
 * All functionality of lsst::afw::image::Exposure is exposed.
 */
template<typename ImageT, typename MaskT=lsst::afw::image::MaskPixel,
         typename VarianceT=lsst::afw::image::VariancePixel>
class CharacterizedExposure : public lsst::afw::image::Exposure<ImageT,MaskT,VarianceT> {
public:
    typedef lsst::afw::image::MaskedImage<ImageT, MaskT, VarianceT> MaskedImageT;
    typedef lsst::afw::image::Exposure<ImageT, MaskT, VarianceT> ExposureT;
    typedef lsst::meas::algorithms::PSF PSF;

    typedef boost::shared_ptr<CharacterizedExposure> Ptr;
    typedef boost::shared_ptr<CharacterizedExposure const> ConstPtr;

    /**
     * Construct a new exposure with the given dimensions with the given wcs,
     * and psf
     */
    explicit CharacterizedExposure(
        int const cols=0, 
        int const rows=0, 
        Wcs const& wcs=Wcs(),
        PSF::Ptr const & psf=PSF::Ptr()
    ) : ExposureT(cols,rows,wcs), _psf(psf) {}

    /**
     * Construct from an existing lsst::afw::image::MaskedImage with the given
     * wcs, and psf
     */
    explicit CharacterizedExposure(
        MaskedImageT & maskedImage, 
        Wcs const& wcs=Wcs(), 
        PSF::Ptr const & psf=PSF::Ptr()
    ) : ExposureT(maskedImage,wcs), _psf(psf) {}

    /**
     * Construct from an existing lsst::afw::image::Exposure by adding a psf
     */
    explicit CharacterizedExposure(
        ExposureT & exposure, 
        PSF::Ptr const & psf=PSF::Ptr()
    ) : ExposureT(exposure), _psf(psf) {}
    
    /**
     * Construct as a subimage of an existing CharacterizedExposure.
     * 
     * @param src master-image to take a sub-image of
     * @param bbox specifies the offset and dimensions of the sub image within
     *      the parent image
     * @param deep specifies whether a deep copy should be performed
     */
    CharacterizedExposure(
        CharacterizedExposure const& src, 
        lsst::afw::image::BBox const& bbox, 
        bool const deep=false
    ) : ExposureT(src,bbox,deep), _psf(src.getPSF()) {}
     
    bool hasPSF() const { return _psf; }
    PSF::Ptr getPSF() const { return _psf; }
    void setPSF(PSF::Ptr const & psf) { _psf = psf; }

private:
    lsst::meas::algorithms::PSF::Ptr _psf;
};

}}}


#endif // !LSST_MEAS_MULTIFIT_CHARACTERIZED_EXPOSURE_H
