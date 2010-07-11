// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
/**
 * @file
 *
 * Declaration of the base Model class
 */
#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include <Eigen/Core>

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>

#include <ndarray_fwd.hpp>

#include "lsst/meas/multifit/core.h"
#include "lsst/afw/geom/Box.h"

namespace lsst{
namespace meas {
namespace multifit {

//forward declarations
class ModelProjection;
class ModelFactory;

/**
 *  Base class for models of astronomical objects.
 *
 *  A Model's parameters are always in "global" (typically celestial) 
 *  coordinates. A projection of the Model to a particular image coordinate 
 *  system (along with convolution by the appropriate PSF and application of 
 *  other observational effects) is represented by an instance of 
 *  ModelProjection, which will generally be subclassed in tandem with Model 
 *  and ModelFactory.  Model and ModelProjection participate in a Observer 
 *  design pattern, with a Model broadcasting changes in its parameters to its
 *  associated ModelProjection objects.
 *
 *  @sa ModelFactory
 *  @sa ModelProjection
 */
class Model : public boost::enable_shared_from_this<Model> {
public:
    typedef boost::shared_ptr<Model> Ptr;
    typedef boost::shared_ptr<Model const> ConstPtr;

    /**
     *  Create a Footprint that would contain a projection of this Model.
     */
    virtual Footprint::Ptr computeProjectionFootprint(
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs
    ) const = 0;

    /**
     *  Create an image-coordinate bounding box that would contain a projection
     *  of this Model.
     */
    virtual lsst::afw::geom::BoxD computeProjectionEnvelope(
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs
    ) const = 0;

    /**
     *  Create an ra/dec bounding ellipse for this Model.
     */
    virtual boost::shared_ptr<lsst::afw::geom::ellipses::Ellipse> computeBoundingEllipse() const = 0;

    /**
     * Immutable access to this Model's linear parameters 
     */
    ParameterVector const & getLinearParameters() const { 
        return *_linearParameters; 
    }

    /**
     * Immutable access to this Model's linear parameters 
     *
     * @return a constant iterator pointing at the beginning of the parameter
     *      vector
     */
    ParameterConstIterator getLinearParameterIter() const { 
        return _linearParameters->data(); 
    }

    /**
     * Number of linear parameters in this Model
     */
    int const getLinearParameterSize() const { 
        return _linearParameters->size(); 
    }

    void setLinearParameters(ParameterConstIterator parameterIter);
    void setLinearParameters(ParameterVector const & parameters);

    /**
     * Immutable access to this Model's nonlinear parameters 
     */
    ParameterVector const & getNonlinearParameters() const { 
        return *_nonlinearParameters; 
    }

    /**
     * Immutable access to this Model's nonlinear parameters 
     *
     * @return a constant iterator pointing at the beginning of the parameter
     *      vector
     */
    ParameterConstIterator getNonlinearParameterIter() const { 
        return _nonlinearParameters->data(); 
    }

    /**
     * Number of nonlinear parameters in this Model.
     */
    int const getNonlinearParameterSize() const { 
        return _nonlinearParameters->size(); 
    }

    void setNonlinearParameters(ParameterConstIterator parameterIter);
    void setNonlinearParameters(ParameterVector const & parameters);

    /**
     *  Create a new Model with the same type and parameter values.
     *
     *  Associated ModelProjection objects are not copied or shared;
     *  the new Model will not have any associated ModelProjections.
     */
    virtual Model::Ptr clone() const = 0;
   
    virtual ~Model() {}

    /** 
     *  Create a ModelProjection object associated with this.
     */
    virtual boost::shared_ptr<ModelProjection> makeProjection(
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint
    ) const = 0;

protected:
    typedef boost::weak_ptr<ModelProjection> ProjectionWeakPtr;
    typedef std::list<ProjectionWeakPtr> ProjectionList;

    /**
     * Initialize the Model and allocate space for the parameter vectors.
     */
    Model(int linearParameterSize, int nonlinearParameterSize)       
      : _linearParameters(), _nonlinearParameters(), _projectionList()
    {
        //Eigen does not handle zero-size matrix creation gracefully.
        if(linearParameterSize > 0)
            _linearParameters.reset(new ParameterVector(linearParameterSize));
        else _linearParameters.reset(new ParameterVector());

        if(nonlinearParameterSize > 0)
            _nonlinearParameters.reset(new ParameterVector(nonlinearParameterSize));
        else _nonlinearParameters.reset(new ParameterVector());

    }

    /**
     * Deep-copy the Model.
     *
     * This is a deep copy of the model parameters, projections will not be
     * associated with the new Model
     */ 
    explicit Model(Model const & model) 
      : _linearParameters(new ParameterVector(model.getLinearParameters())),
        _nonlinearParameters(new ParameterVector(model.getNonlinearParameters())),
        _projectionList()
    {}

    explicit Model() :
       _linearParameters(), 
       _nonlinearParameters(), 
       _projectionList() 
    {}

    void _broadcastLinearParameterChange() const;

    void _broadcastNonlinearParameterChange() const;

    /**
     *  Provide additional handling code for setting linear parameter.
     *
     *  This will be called by setLinearParameters, after the parameter 
     *  vector has been updated and before the call to 
     *  _broadcastLinearParameterChange().
     */
    virtual void _handleLinearParameterChange() {}

    /**
     *  Provide additional handling code for setting nonlinear parameters.
     *
     *  This will be called by setNonlinearParameters, after the parameter 
     *  vector has been updated and before the call to 
     *  _broadcastNonlinearParameterChange().
     */
    virtual void _handleNonlinearParameterChange() {}

    /**
     * Add a newly-created projection to the list of listeners.
     */
    void _registerProjection(
        boost::shared_ptr<ModelProjection> const & projection
    ) const;

    boost::shared_ptr<ParameterVector> _linearParameters;
    boost::shared_ptr<ParameterVector> _nonlinearParameters;

private:
    friend class ModelFactory;

    // disable default constructor and assignment

    void operator=(Model const & other) { assert(false); } 

    mutable ProjectionList _projectionList;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_H
