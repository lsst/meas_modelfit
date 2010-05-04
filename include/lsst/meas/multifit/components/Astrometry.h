// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "lsst/afw/image/Utils.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModel;
class ComponentModelProjection;

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
    enum Parameters {RA=0, DEC, SIZE};

    typedef boost::shared_ptr<Astrometry> Ptr;
    typedef boost::shared_ptr<Astrometry const> ConstPtr;
    typedef Eigen::Matrix<Pixel, SIZE, Eigen::Dynamic> DerivativeMatrix;

    virtual ~Astrometry() {}

    Astrometry() 
      : _parameters(new ParameterVector(SIZE)), 
        _parameterIter(_parameters->data()) 
    {}

    /** 
     * Copy construct an Astrometry object
     *
     * This is by a shallow copy by default. If deep ==true, a deep copy will be
     * performed
     */
    /// Construct an Astrometry object for use a fixed position
    explicit Astrometry(lsst::afw::geom::Point2D const & position) : 
        _parameters(new ParameterVector(position.asVector())),
        _parameterIter(_parameters->data())
    {}
    

    explicit Astrometry(Astrometry const & other) :
        _parameters(new ParameterVector(other.computePosition().asVector())),
        _parameterIter(_parameters->data())
    {}
    
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
            _parameterIter[RA], _parameterIter[DEC]
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



protected:
    friend class multifit::ComponentModel;
    friend class multifit::ComponentModelProjection;

    /// Return the number of astrometric parameters.
    virtual int const getParameterSize() const { return SIZE; }


    /**
     * Handle a change in the parameters.
     */
    virtual void _handleParameterChange() {}

    /**
     * Construct a new Astrometry object, using this as a type factory
     * 
     * @sa Astrometry(boost::shared_ptr<ParameterVector> const&, size_t const &, bool const &)
     */
    virtual Ptr create (
        boost::shared_ptr<ParameterVector> const & parameters, 
        size_t const & start = 0
    ) const {
        return Ptr(new Astrometry(parameters, start));    
    }

    ParameterConstIterator begin() const {return _parameterIter;}
    ParameterConstIterator end() const {return _parameterIter + getParameterSize();}


private:
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
        boost::shared_ptr<ParameterVector const> const & parameters, 
        size_t const & start=0        
    ) : _parameters(parameters), 
        _parameterIter(_parameters->data() + start)
    {
        if(!parameters) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Uninitialized input parameter vector"
            );
        }

        if(static_cast<int>(start + SIZE) >  parameters->size()) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException, 
                "Input parameter vector must be at least of length start + getParameterSize()"
            );
        }
    }
    boost::shared_ptr<ParameterVector const> _parameters;
    ParameterConstIterator _parameterIter;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_ASTROMETRY_H

