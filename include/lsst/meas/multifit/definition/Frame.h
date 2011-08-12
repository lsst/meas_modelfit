// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_Frame
#define LSST_MEAS_MULTIFIT_DEFINITION_Frame

#include "lsst/ndarray.h"
#include "lsst/meas/multifit/constants.h"
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include "boost/serialization/binary_object.hpp"

namespace boost {
namespace serialization {
    class access;
}}

namespace lsst { namespace meas { namespace multifit {

namespace detail {

/// @brief Utility (nonpolymorphic) base class for definition::Frame and grid::Frame.
class FrameBase {
public:

    ID const id;

#ifndef SWIG

    /// @brief Return the Frame's filter ID.
    FilterId const getFilterId() const { return _filterId; }

    /// @brief Return the Frame's Wcs.
    Wcs::ConstPtr const getWcs() const { return _wcs; }

    /// @brief Return the Frame's Psf.
    Psf::ConstPtr const getPsf() const { return _psf; }

    /// @brief Return the footprint that defines the pixel region of interest.
    CONST_PTR(Footprint) const getFootprint() const { return _footprint; }

    /// @brief Return the 1-d data array (with size matching the footprint area).
    lsst::ndarray::Array<Pixel const,1,1> const getData() const { return _data; }

    /// @brief Return the 1-d pixel weight array (with size matching the footprint area).
    lsst::ndarray::Array<Pixel const,1,1> const getWeights() const { return _weights; }

#endif

protected:
    
    FrameBase(
        ID id_,
        Footprint::Ptr const & footprint,
        lsst::ndarray::Array<Pixel,1,1> const & data,
        lsst::ndarray::Array<Pixel,1,1> const & weights = lsst::ndarray::Array<Pixel,1,1>()
    ) : id(id_), _footprint(footprint), _data(data), _weights(weights) {}

    FrameBase(FrameBase const & other, bool deep=false);
    FrameBase(ID id_) : id(id_) {};
    FilterId _filterId;
    Wcs::Ptr _wcs;
    Psf::Ptr _psf;
    Footprint::Ptr _footprint;
    lsst::ndarray::Array<Pixel,1,1> _data;
    lsst::ndarray::Array<Pixel,1,1> _weights;
};

} // namespace detail


namespace definition {

/**
 *  @brief A customized exposure-like class for multifit definitions.
 *
 *  Like those of definition::ObjectComponent, accessors of Frame return by non-const reference, reflecting
 *  the fact that they can be set freely and will be validated only when a Grid is constructed.
 *  from the Definition.  We have used accessors rather than public data members so the interface 
 *  is similar to that of grid::Frame, which behaves like a const version of definition::Frame.
 */
class Frame : public detail::FrameBase {
public:
    
    Frame(
        ID id_, 
        Footprint::Ptr const & footprint,
        lsst::ndarray::Array<Pixel,1,1> const & data,
        lsst::ndarray::Array<Pixel,1,1> const & weights = lsst::ndarray::Array<Pixel,1,1>()
    ) : detail::FrameBase(id_, footprint, data, weights) {}

    Frame(
        ID id_,
        FilterId filterId,
        Wcs::Ptr const & wcs,
        Psf::Ptr const & psf,
        Footprint::Ptr const & footprint,
        lsst::ndarray::Array<Pixel,1,1> const & data,
        lsst::ndarray::Array<Pixel,1,1> const & weights = lsst::ndarray::Array<Pixel,1,1>()
    ) : detail::FrameBase(id_, footprint, data, weights) {
        _filterId = filterId; _wcs = wcs; _psf = psf;
    }

    Frame(Frame const & other) : detail::FrameBase(other, false) {}

    /**
     *  @brief Construct a frame from an exposure and footprint.
     *
     *  Masked pixels will be removed from the footprint using
     *  afw::detection::footprintAndMask with the given bitmask.
     *  The default value is to remove any masked pixels.
     */
    template <typename PixelT>
    static Frame make(
        ID const id,
        lsst::afw::image::Exposure<PixelT> const & exposure,
        Footprint::Ptr const & footprint,
        lsst::afw::image::MaskPixel const bitmask=~0x0
    );
    
#ifndef SWIG

    /// @brief Return and/or set the Frame's filter ID.
    FilterId & getFilterId() { return _filterId; }
    using detail::FrameBase::getFilterId;

    /// @brief Return and/or set the Frame's Wcs.
    Wcs::Ptr & getWcs() { return _wcs; }
    using detail::FrameBase::getWcs;

    /// @brief Return and/or set the Frame's Psf.
    Psf::Ptr & getPsf() { return _psf; }
    using detail::FrameBase::getPsf;

    /// @brief Return and/or set the footprint that defines the pixel region of interest.
    Footprint::Ptr & getFootprint() { return _footprint; }
    using detail::FrameBase::getFootprint;

    /// @brief Return and/or set the 1-d data array (with size matching the footprint area).
    lsst::ndarray::Array<Pixel,1,1> & getData() { return _data; }
    using detail::FrameBase::getData;

    /// @brief Return and/or set the 1-d pixel weight array (with size matching the footprint area).
    lsst::ndarray::Array<Pixel,1,1> & getWeights() { return _weights; }
    using detail::FrameBase::getWeights;

#endif

    //@{
    /**
     *  @brief Attribute setters.
     *
     *  These are unnecessary in C++ because the non-const getters return by reference,
     *  but are included here because they're the only way to set thing in Python.
     */
    void setFilterId(FilterId filterId) { _filterId = filterId; }
    void setWcs(Wcs::Ptr const & wcs) { _wcs = wcs; }
    void setPsf(Psf::Ptr const & psf) { _psf = psf; }
    void setFootprint(Footprint::Ptr const & footprint) { _footprint = footprint; }
    void setData(lsst::ndarray::Array<Pixel,1,1> const & data) { _data = data; }
    void setWeights(lsst::ndarray::Array<Pixel,1,1> const & weights) { _weights = weights; }
    //@}

private:
    friend class grid::Initializer;
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, unsigned int const version) {}
    explicit Frame(detail::FrameBase const & other) : detail::FrameBase(other, false) {}
};



#ifndef SWIG
std::ostream & operator<<(std::ostream & os, Frame const & frame);
#endif

}}}} // namespace lsst::meas::multifit::definition

namespace boost { namespace serialization {    
template <class Archive>
inline void save_construct_data(
    Archive & ar, 
    const lsst::meas::multifit::definition::Frame * frame, 
    unsigned int const version
) {
    ar << make_nvp("id", frame->id);
    ar << make_nvp("filterId", frame->getFilterId());
    ar << make_nvp("wcs", frame->getWcs());
    ar << make_nvp("psf", frame->getPsf());
    lsst::meas::multifit::Footprint::ConstPtr fp = frame->getFootprint();
    ar << make_nvp("footprint", fp);

    int size = fp->getArea();
    ar << make_nvp("data", make_array(frame->getData().getData(), size));
    ar << make_nvp("weights", make_array(frame->getWeights().getData(), size));
}

template <class Archive>
inline void load_construct_data(    
    Archive & ar, 
    lsst::meas::multifit::definition::Frame * frame, 
    unsigned int const version
) {
    lsst::meas::multifit::ID id;
    lsst::meas::multifit::FilterId filterId;
    lsst::meas::multifit::Psf::Ptr psf;
    lsst::meas::multifit::Wcs::Ptr wcs;
    lsst::meas::multifit::Footprint::Ptr fp;
    lsst::ndarray::Array<lsst::meas::multifit::Pixel, 1, 1> data, weights;

    ar >> make_nvp("id", id);
    ar >> make_nvp("filterId", filterId);
    ar >> make_nvp("wcs", wcs);
    ar >> make_nvp("psf", psf);
    ar >> make_nvp("footprint", fp);

    int size = fp->getArea();
    data = lsst::ndarray::allocate(size);
    weights = lsst::ndarray::allocate(size);
    ar >> make_nvp("data", make_array(data.getData(), size));
    ar >> make_nvp("weights", make_array(weights.getData(), size));


    ::new(frame) lsst::meas::multifit::definition::Frame(
        id, filterId, wcs, psf, fp, data, weights
    );
}

}} //namespace boost::serialization
#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Frame
