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
 * Includes declaration of ModelEvaluator
 */
#ifndef LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H
#define LSST_MEAS_MULTIFIT_MODEL_EVALUATOR_H

#include <vector>
#include <iostream>
#include "Eigen/Core"
#include "Eigen/LU"

#include "ndarray.hpp"

#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/afw/geom/AffineTransform.h"

namespace lsst {
namespace meas {
namespace multifit{

/**
 * Manage projection of a model to a set of Exposures
 *
 * The ModelEvaluator abstracts the notion of multifit by composing
 * grand-matrices which are populated, part by part, by
 * ModelProjection objects created for each Exposure. 
 *
 */
class ModelEvaluator : private boost::noncopyable {
public:
    typedef std::list<ModelProjection::Ptr> ProjectionList;
    typedef boost::shared_ptr<ModelEvaluator> Ptr;
    typedef boost::shared_ptr<const ModelEvaluator> ConstPtr;

    /**
     * Construct a ModelEvaluator
     *
     * @param model Model to manage
     * @param nMinPix minimum number of pixels an exposure must contribute to 
     *      the model's projected footprint to be used
     */
    explicit ModelEvaluator(
        Model::ConstPtr const & model, 
        lsst::afw::geom::AffineTransform const &pixelToSky,
        int const nMinPix=0) 
      : _validProducts(0),
        _model(model->clone()),
        _pixelToSky(pixelToSky)
    {
        setMinPixels(nMinPix);
    }

    template <typename ExposureT>
    void setExposures(
        std::list<ExposureT> const & exposureList
    );

    template <typename MaskedImageT>
    void setData(
        std::list<MaskedImageT> const & imageList,
        std::list<lsst::afw::detection::Psf::ConstPtr> const & psfList,
        std::list<lsst::afw::geom::AffineTransform> const & skyToPixelList
    );

    template <typename MaskedImageT>
    void setData(
        MaskedImageT const  & image, 
        CONST_PTR(lsst::afw::detection::Psf) const & psf,
        lsst::afw::geom::AffineTransform const &skyToPixel
    );
#ifndef SWIG
    /**
     * Vector of image data from all contributing pixels from all the exposures
     *
     * The data vector is composed from the concactenation of each exposure's
     * footprint-compressed image data. The resulting data vector is an
     * abstraction of all the pixels this model is evaluated on.
     *
     * @sa getVarianceVector
     */
    ndarray::Array<Pixel const, 1, 1> getDataVector() const {
        return _dataVector;
    }
#endif
    /**
     * Compute the sigma for each contributing pixel from all exposure's
     *
     * The sigma vector is the component wise sqaure root of the variance.
     *
     * @sa getVarianceVector
     */
    Eigen::Matrix<Pixel, Eigen::Dynamic, 1> const & getSigmaVector() const {
        return _sigma;
    }

    Eigen::Matrix<Pixel, Eigen::Dynamic, 1> const getWeightedData() const {
        VectorMap map(_dataVector.getData(), getNPixels(), 1);
        return map.cwise()/_sigma;
    }

    /**
     * @name Model Product Computers
     *     
     * Each of these functions compute a footprint-compressed product as a
     * column-major matrix with number of rows corresponding to the
     * footprint-mapped pixel index, and the number of columns 
     * corresponding to the parameter index.
     */
    //@{
    Eigen::Matrix<Pixel, Eigen::Dynamic, 1> const & computeModelImage();
    Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> const & computeLinearParameterDerivative();
    Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> const & computeNonlinearParameterDerivative();
    //@}

    /**
     * Pixel threshold used to discriminate exposures
     *
     * Only exposures on which the Model's footprint covers more than this
     * threshold will be used to evaluate the model
     */
    int const & getMinPixels() const {return _nMinPix;}

    /**
     * Pixel threshold used to discriminate exposures
     *
     * Only exposures on which the Model's footprint covers more than this
     * threshold will be used to evaluate the model.
     *
     * To disable discrimination, call with nMinPIx <= 0
     *
     * @param nMinPix new threshold to to use to discard exposures with too few
     *      contributing pixels. if nMinPix <= 0, no exposure's will be
     *      discarded
     */
    void setMinPixels(int const nMinPix) {
        _nMinPix = (nMinPix < 0)? 0 : nMinPix;
    }

    /**
     * @name Model Accessors
     *
     * Vet access to the model by exposing a subset of Model's functionality
     * directly
     */
    //@{
    /**
     * Number of linear parameters in the model
     */
    int const getLinearParameterSize() const {
        return _model->getLinearParameterSize();
    }
    /**
     * Number of nonlinear parameters in the model
     */
    int const getNonlinearParameterSize() const {
        return _model->getNonlinearParameterSize();
    }
    /**
     * Immutable access to the model's linear parameters
     */
    ParameterVector const & getLinearParameters() const {
        return _model->getLinearParameters();
    }
    /**
     * Immutable access to the model's nonlinear parameters
     */
    ParameterVector const & getNonlinearParameters() const {
        return _model->getNonlinearParameters();
    }

    /**
     * Set the model's linear parameters
     */
    void setLinearParameters(ParameterVector const & parameters) {
        setLinearParameters(parameters.data());
    }

    /**
     * Set the model's linear parameters
     */
    void setLinearParameters(
        ParameterConstIterator const & parameterIterator
    ) {    
        _model->setLinearParameters(parameterIterator);
        _validProducts &= (~MODEL_IMAGE);
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }

    /**
     * Set the model's nonlinear parameters
     */
    void setNonlinearParameters(ParameterVector const & parameters) {
        setNonlinearParameters(parameters.data());
    }

    /**
     * Set the model's nonlinear parameters
     */
    void setNonlinearParameters(
        ParameterConstIterator const & parameterIterator
    ) {
        _model->setNonlinearParameters(parameterIterator);
        _validProducts &= (~MODEL_IMAGE);
        _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }
    //@}  
    
    /**
     * Immutable access to the model this evaluator is managing
     */
    Model::ConstPtr getModel() const {return _model;}
    /**
     * number of exposures used
     */
    int const getNProjections() const {return _projectionList.size();}
    /**
     * Number of total pixels accross all exposures used
     */
    int const getNPixels() const {return _dataVector.getSize<0>();}

    /**
     * List the ModelProjection objects
     *
     * The number of projections may not match the number of input exposures
     * as the evluator discards any exposures with fewer than nMinPix
     * pixels covered by the model's projected footprint.
     */
    ProjectionList const & getProjectionList() const {
        return _projectionList;
    }

private:   

    typedef ProjectionList::iterator ProjectionIterator;
    typedef lsst::afw::detection::Footprint::SpanList SpanList;


    enum ProductFlag {
        MODEL_IMAGE = 1<<0,
        LINEAR_PARAMETER_DERIVATIVE = 1<<1,
        NONLINEAR_PARAMETER_DERIVATIVE = 1<<2,
    };

    int _nMinPix;
    int _validProducts;
    Model::Ptr _model;    
    ProjectionList _projectionList;
    lsst::afw::geom::AffineTransform _pixelToSky;
    
    ndarray::Array<Pixel, 1, 1> _dataVector;
    ndarray::Array<Pixel, 1, 1> _varianceVector;    
    ndarray::Array<Pixel, 1, 1> _modelImageBuffer;
    ndarray::Array<Pixel, 2, 2> _linearDerivativeBuffer;
    ndarray::Array<Pixel, 2, 2> _nonlinearDerivativeBuffer;

    Eigen::Matrix<Pixel, Eigen::Dynamic, 1> _modelImage, _sigma;
    Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> _linearDerivative, _nonlinearDerivative;


};

}}} //end namespace lsst::meas::multifit

#endif
