#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "lsst/afw/image/Utils.h"

#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModel;

namespace components {

/**
 *  \brief A concrete base class for component objects that decribe the position of 
 *  a model.
 *
 *  The base class has only two parameters, the RA and DEC of the model.  Subclasses
 *  will generally treat this points as a "reference point" and add additional parameters
 *  that compute the output point relative to the reference point.
 */
class Astrometry {
public:
    typedef boost::shared_ptr<Astrometry> Ptr;
    typedef boost::shared_ptr<Astrometry const> ConstPtr;

    typedef Eigen::Matrix<Pixel,2,Eigen::Dynamic> DerivativeMatrix;

    virtual ~Astrometry() {}

    /// \brief Default-construct an Astrometry object for use as a template.
    Astrometry() : _astrometryParameterIter(NULL) {}

    /// \brief Return the number of astrometric parameters.
    virtual int const getAstrometryParameterSize() const { return 2; }

    /// \brief Return the reference point (ra,dec).
    lsst::afw::geom::Point2D getRefPoint() const { 
        return lsst::afw::geom::Point2D::makeXY(_astrometryParameterIter[0],_astrometryParameterIter[1]);
    }

    /**
     *  Return the computed point (ra,dec).
     */
    virtual lsst::afw::geom::Point2D apply() const { return getRefPoint(); }

    /**
     *  Compute the derivative with respect to the astrometric parameters.
     */
    virtual DerivativeMatrix const & differentiate() const {
        static DerivativeMatrix i = DerivativeMatrix::Identity(2,2);
        return i;
    }

protected:

    friend class multifit::ComponentModel;

    /**
     *  \brief Construct a new Astrometry object using this as a template.
     *
     *  Typically only used by ComponentModel.  The passed iterator must
     *  remain valid for the full lifetime of the Astrometry object.
     *
     *  \param astrometryParameterIter pointer to the first Astrometry-specific
     *   parameter in the owning ComponentModel's nonlinear parameter vector 
     *   
     */
    virtual Astrometry::Ptr create(
        ParameterConstIterator astrometryParameterIter 
    ) const {
        return Astrometry::Ptr(new Astrometry(astrometryParameterIter));
    }

    virtual void _handleAstrometryParameterChange() {}

private:

    /**
     *  \brief Construct an Astrometry object for use inside a ComponentModel.
     *
     *  \sa Astrometry::create()
     */
    explicit Astrometry(
        ParameterConstIterator astrometryParameterIter
    ) : _astrometryParameterIter(astrometryParameterIter) {}

    void operator=(Astrometry const & other) { assert(false); } // Assignment disabled.

    ParameterConstIterator _astrometryParameterIter; 
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H

