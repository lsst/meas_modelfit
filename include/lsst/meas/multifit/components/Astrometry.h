// -*- lsst-c++ -*-
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
 *  Concrete base class for component objects that decribe the position of 
 *  a model.
 *
 *  The base class has only two parameters, the RA and DEC of the model.  
 *  Subclasses will generally treat this points as a "reference point" and add 
 *  additional parameters that compute a position relative to the reference 
 *  point.
 */
class Astrometry {
public:
    typedef boost::shared_ptr<Astrometry> Ptr;
    typedef boost::shared_ptr<Astrometry const> ConstPtr;
    typedef Eigen::Matrix<Pixel, 2, Eigen::Dynamic> DerivativeMatrix;

    virtual ~Astrometry() {}

    /// Default-construct an Astrometry object for use as a template.
    Astrometry() : _astrometryParameterIter(NULL) {}

    /// Return the number of astrometric parameters.
    virtual int const getParameterSize() const { return 2; }

    /**
     * Return the reference point (ra,dec).
     *
     * For static sources, the computed position is equivalent to the reference
     * point
     *
     * @sa computePosition
     */
    lsst::afw::geom::Point2D getReferencePoint() const { 
        return lsst::afw::geom::Point2D::make(
            _astrometryParameterIter[0], _astrometryParameterIter[1]
        );
    }

    /**
     *  Return the computed point (ra,dec).
     *
     *  For static sources, the computed position is equivalent to the reference
     *  point
     *
     *  @sa getReferencePoint
     */
    virtual lsst::afw::geom::Point2D computePosition() const { 
        return getReferencePoint(); 
    }

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
     *  Construct a new Astrometry object using this as a template.
     *
     *  Only used by ComponentModel. The passed iterator must
     *  remain valid for the full lifetime of the Astrometry object.
     *
     *  @param astrometryParameterIter pointer to the first Astrometry-specific
     *   parameter in the owning ComponentModel's nonlinear parameter vector 
     *   
     */
    virtual Astrometry::Ptr create(
        ParameterConstIterator astrometryParameterIter 
    ) const {
        return Astrometry::Ptr(new Astrometry(astrometryParameterIter));
    }

    /**
     * Handle a change in the parameters.
     */
    virtual void _handleParameterChange() {}

private:

    /**
     *  Construct an Astrometry object for use inside a ComponentModel.
     *
     *  @sa Astrometry::create()
     */
    explicit Astrometry(
        ParameterConstIterator astrometryParameterIter
    ) : _astrometryParameterIter(astrometryParameterIter) {}

    //disable asignment
    void operator=(Astrometry const & other) { assert(false); }

    ParameterConstIterator _astrometryParameterIter; 
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H

