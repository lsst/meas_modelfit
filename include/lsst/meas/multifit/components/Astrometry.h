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

    /** 
     * Copy construct an Astrometry object
     *
     * This is by a shallow copy by default. If deep ==true, a deep copy will be
     * performed
     */
    explicit Astrometry(
        Astrometry const & other, 
        bool const & deep=false
    ) :  _parameters(other._parameters),
        _astrometryParameterIter(other._astrometryParameterIter)
    {
        if(!_parameters) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::NullPointerException,
                "Uninitialized parameter vector"
            );
        }
        if(deep) {
            Astrometry * tmp = new Astrometry(computePosition());
            swap(this, tmp);
            delete tmp;
        }
    }

    /// Construct an Astrometry object for use a fixed position
    explicit Astrometry(lsst::afw::geom::Point2D const & position) : 
        _parameters(new ParameterVector(position.asVector())),
        _astrometryParameterIter(_parameters->begin())
    {}
    
    /** 
     * Construct an Astrometry object to use to interpret a range of a 
     * ParameterVector
     *
     * @param parameters must have length at least start+getParameterSize()
     * @param start vector index of the first astrometry relevant parameter
     *
     * @throws lsst::pex::exceptions::InvalidParameterException if parameters is
     * not at least of length start+getParameterSize()
     */
    explicit Astrometry(
        boost::shared_ptr<ParameterVector> const & parameters, 
        size_t const & start=0,
        bool deep = false
    ) : _parameters(parameters), 
        _astrometryParameterIter(_parameters->begin() + start)
    {
        if(!_parameters) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::NullPointerException,
                "Uninitialized input parameter vector"
            );
        }

        if(start + SIZE >  parameters->size()) {
            throw LSST_EXCEPT(
                lsst::pex::exception::InvalidParameteException, 
                "Input parameter vector must be at least of length start + getParameterSize()"
            );
        }
        if(deep)
    }
    
    /// Return the number of astrometric parameters.
    virtual int const getParameterSize() const { return SIZE; }

    /**
     * Return the reference point (ra,dec).
     *
     * For static sources, the computed position is equivalent to the reference
     * point
     *
     * @throw lsst::pex::exceptions::NullPointerException if not initialized
     * @sa computePosition
     */
    virtual lsst::afw::geom::Point2D getReferencePoint() const { 
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
        static DerivativeMatrix i = DerivativeMatrix::Identity(SIZE,SIZE);
        return i;
    }

    /**
     * Handle a change in the parameters.
     */
    virtual void handleParameterChange() {}

    /**
     *  Construct a new Astrometry object using this as a template.
     * 
     * @sa Astrometry(boost::shared_ptr<ParameterVector const> const&, size_t const &, bool const &)
     */
    virtual Ptr create(
        boost::shared_ptr<ParameterVector const> const & parameters, 
        size_t const & start = 0
    ) {
        return boost::make_shared<Astrometry>(parameters, start, false);
    }
    
    virtual Ptr constrain() const;
protected:
    void deepCopy(Astrometry const & other) {
        _parameters.reset(new ParameterVector(other.computePosition().asVector()));
        _astrometryParameterIter=_parameters.begin();
    }

    static const int SIZE = 2;
    Astrometry() : _parameters(NULL), _astrometryParameterIter(NULL) {}

    boost::shared_ptr<ParameterVector const> _parameters;
    ParameterConstIterator _astrometryParameterIter;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H

