// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Model

#include "boost/pointer_cast.hpp"
#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"

namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;

class DoNothingProjection : public multifit::ModelProjection {
public:
    typedef boost::shared_ptr<DoNothingProjection> Ptr;

    DoNothingProjection(
        multifit::Model::ConstPtr const & model, 
        multifit::WcsConstPtr const & wcs, 
        multifit::FootprintConstPtr const & footprint
    ) : ModelProjection(model, wcs, footprint), _linearChange(0), _nonlinearChange(0) {}

    virtual int const getWcsParameterSize() const {return 0;}
    virtual int const getPsfParameterSize() const {return 0;}
   
    virtual bool hasWcsParameterDerivative() const {return false;}
    virtual bool hasPsfParameterDerivative() const {return false;}

    int getLinearChange() const {return _linearChange;}
    int getNonlinearChange() const {return _nonlinearChange;}
protected:
    virtual void _computeModelImage(ndarray::Array<multifit::Pixel, 1, 1> const & vector) {}
    virtual void _computeLinearParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    virtual void _computeNonlinearParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    virtual void _computePsfParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    virtual void _computeWcsParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    
    virtual void _handleLinearParameterChange() {++_linearChange;}
    virtual void _handleNonlinearParameterChange() {++_nonlinearChange;}

    int _linearChange, _nonlinearChange;
};

class DoNothingModel : public multifit::Model {
public:
    //initialized a model with 1 linear, and 1 nonlinear parameter
    typedef boost::shared_ptr<DoNothingModel> Ptr;

    static DoNothingModel::Ptr makeModel() {
        return DoNothingModel::Ptr(new DoNothingModel());
    }
    
    virtual multifit::Footprint::Ptr computeProjectionFootprint(
        multifit::PsfConstPtr const & psf,
        multifit::WcsConstPtr const & wcs
    ) const {
        return boost::make_shared<multifit::Footprint>(
            lsst::afw::image::BBox(
                lsst::afw::image::PointI(0, 0),
                10, 15
            )
        );
    }
    virtual lsst::afw::geom::BoxD computeProjectionEnvelope(
        multifit::PsfConstPtr const & psf,
        multifit::WcsConstPtr const & wcs
    ) const {
        return lsst::afw::geom::BoxD(
            lsst::afw::geom::PointD::make(0, 0),
            lsst::afw::geom::ExtentD::make(10, 15)
        );
    }
    virtual lsst::afw::geom::ellipses::Ellipse::Ptr computeBoundingEllipse() const {
        lsst::afw::geom::ellipses::Axes core(8, 6, 0);
        return boost::make_shared<lsst::afw::geom::ellipses::AxesEllipse>(
            core,
            lsst::afw::geom::PointD::make(5, 7.5)
        );
    }

    virtual multifit::ModelProjection::Ptr makeProjection(
        multifit::PsfConstPtr const & psf,
        multifit::WcsConstPtr const &wcs,
        multifit::FootprintConstPtr const & footprint
    ) const {
        multifit::ModelProjection::Ptr projection(
            new DoNothingProjection(shared_from_this(),wcs, footprint)
        );
        _registerProjection(projection);
        return projection;
    }

    multifit::Model::Ptr clone() const {
        return makeModel();
    }


protected:

    //initialize a model with 1 linear, 1 nonlinear parameter
    DoNothingModel() : multifit::Model(1, 1) {}
};

    
BOOST_AUTO_TEST_CASE(ModelBasic) {

    multifit::Model::Ptr model = DoNothingModel::makeModel();
    BOOST_CHECK_EQUAL(model->getLinearParameterSize(), 1);
    BOOST_CHECK_EQUAL(model->getNonlinearParameterSize(), 1);
    BOOST_CHECK(model->getLinearParameterIter() != 0);
    BOOST_CHECK(model->getNonlinearParameterIter() != 0);


    multifit::ParameterVector linear(1), nonlinear(1);
    linear[0] = 5;
    nonlinear[0] = 6;
    BOOST_CHECK_NO_THROW(model->setLinearParameters(linear.data()));
    BOOST_CHECK_NO_THROW(model->setNonlinearParameters(nonlinear.data()));
    BOOST_CHECK_EQUAL(linear, model->getLinearParameterVector());
    BOOST_CHECK_EQUAL(nonlinear, model->getNonlinearParameterVector());
    
    multifit::WcsConstPtr wcs;
    multifit::PsfConstPtr psf;
    multifit::FootprintConstPtr footprint = model->computeProjectionFootprint(psf, wcs);
    multifit::ModelProjection::Ptr projection = model->makeProjection(psf, wcs, footprint);
    BOOST_CHECK_EQUAL(projection->getModel(), model);
    BOOST_CHECK_EQUAL(projection->getPsfParameterSize(), 0);
    BOOST_CHECK_EQUAL(projection->getWcsParameterSize(), 0);

    DoNothingProjection::Ptr doNothingProjection = boost::static_pointer_cast<DoNothingProjection>(projection);
    BOOST_CHECK_NO_THROW(model->setLinearParameters(linear.data()));
    BOOST_CHECK_NO_THROW(model->setNonlinearParameters(nonlinear.data()));
    BOOST_CHECK_EQUAL(doNothingProjection->getLinearChange(), 1);
    BOOST_CHECK_EQUAL(doNothingProjection->getNonlinearChange(), 1);

    BOOST_CHECK_NO_THROW(model->setLinearParameters(linear.data()));
    BOOST_CHECK_NO_THROW(model->setNonlinearParameters(nonlinear.data()));
    BOOST_CHECK_EQUAL(doNothingProjection->getLinearChange(), 2);
    BOOST_CHECK_EQUAL(doNothingProjection->getNonlinearChange(), 2);
}
