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
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/ComponentModelFactory.h"

namespace multifit = lsst::meas::multifit;
namespace components = lsst::meas::multifit::components;
namespace geom = lsst::afw::geom;

class DoNothingMorphologyProjection : public components::MorphologyProjection {
public:
    DoNothingMorphologyProjection(
        components::Morphology::ConstPtr morphology, 
        geom::ExtentI const & kernelDimensions,
        geom::AffineTransform::ConstPtr const & transform
    ) : components::MorphologyProjection(morphology, kernelDimensions, transform) {}

    virtual ParameterJacobianMatrix const & computeProjectedParameterJacobian() const {
        static const ParameterJacobianMatrix m;
        return m;
    }

    virtual TransformJacobianMatrix const & computeTransformParameterJacobian() const {
        static const TransformJacobianMatrix m;
        return m;
    }
};
class DoNothingMorphology : public components::Morphology {
public:
    typedef boost::shared_ptr<DoNothingMorphology> Ptr;
    virtual int const getMinLinearParameterSize() const {return 1;}
    virtual int const getMaxLinearParameterSize() const {return 1;}
    virtual int const getMorphologyParameterSize() const {return 0;}

    static DoNothingMorphology::Ptr createTemplate() {
        return DoNothingMorphology::Ptr(new DoNothingMorphology());
    }

    virtual geom::ellipses::Core::Ptr computeBoundingEllipseCore() const {
        return boost::make_shared<geom::ellipses::Axes>(8, 6, 0);
    }

    virtual components::MorphologyProjection::Ptr makeProjection(
        geom::Extent2I const & kernelDimensions, 
        geom::AffineTransform::ConstPtr const & transform
    ) const {
        return boost::make_shared<DoNothingMorphologyProjection>(
            shared_from_this(), kernelDimensions, transform
        );
    }

    virtual components::Morphology::Ptr create(
        boost::shared_ptr<multifit::ParameterVector const> const & linearParameterVector,
        multifit::ParameterConstIterator morphologyParameterIter
    ) const {
        return components::Morphology::Ptr(
            new DoNothingMorphology(linearParameterVector, morphologyParameterIter)
        );
    }

protected:
    DoNothingMorphology(
        boost::shared_ptr<multifit::ParameterVector const> const & linearParameterVector,
        multifit::ParameterConstIterator morphologyParameterIter
    ) : components::Morphology(linearParameterVector, morphologyParameterIter) {}


    DoNothingMorphology() : components::Morphology() {}
};

class DoNothingProjection : public multifit::ComponentModelProjection {
public:
    typedef boost::shared_ptr<DoNothingProjection> Ptr;

    DoNothingProjection(
        multifit::ComponentModel::ConstPtr const & model, 
        multifit::PsfConstPtr const & psf,
        multifit::WcsConstPtr const & wcs, 
        multifit::FootprintConstPtr const & footprint
    ) : multifit::ComponentModelProjection(model, psf, wcs, footprint), 
        _linearChange(0), 
        _nonlinearChange(0) 
    {}

    virtual int const getPsfParameterSize() const {return 0;}
    
    int getLinearChange() const {return _linearChange;}
    int getNonlinearChange() const {return _nonlinearChange;}
protected:

    //from ComponentModelProjection
    virtual void _computeTranslationDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    virtual void _computeProjectedParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    
    //from ModelProjection
    virtual void _computeLinearParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    virtual void _computePsfParameterDerivative(ndarray::Array<multifit::Pixel, 2, 1> const & matrix) {}
    
    virtual void _handleLinearParameterChange() {++_linearChange;}
    virtual void _handleNonlinearParameterChange() {++_nonlinearChange;}

    int _linearChange, _nonlinearChange;
};

BOOST_AUTO_TEST_CASE(ModelBasic) {
    multifit::ModelFactory::ConstPtr factory = multifit::ComponentModelFactory::create(
        boost::make_shared<components::Astrometry>(),
        DoNothingMorphology::createTemplate()
    );
    multifit::ParameterVector linear(1), nonlinear(2);
    multifit::Model::Ptr model = factory->makeModel(1, linear.data(), nonlinear.data());
    BOOST_CHECK_EQUAL(model->getLinearParameterSize(), 1);
    BOOST_CHECK_EQUAL(model->getNonlinearParameterSize(), 2);
    BOOST_CHECK(model->getLinearParameterIter() != 0);
    BOOST_CHECK(model->getNonlinearParameterIter() != 0);


    linear[0] = 5;
    nonlinear[0] = 6;
    nonlinear[1] = 7;
    BOOST_CHECK_NO_THROW(model->setLinearParameters(linear.data()));
    BOOST_CHECK_NO_THROW(model->setNonlinearParameters(nonlinear.data()));
    BOOST_CHECK_EQUAL(linear, model->getLinearParameterVector());
    BOOST_CHECK_EQUAL(nonlinear, model->getNonlinearParameterVector());
    
    multifit::WcsConstPtr wcs = boost::make_shared<lsst::afw::image::Wcs>();
    multifit::PsfConstPtr psf =lsst::meas::algorithms::createPSF("DoubleGaussian", 7, 7, 0.33);
    multifit::FootprintConstPtr footprint = model->computeProjectionFootprint(psf, wcs);
    multifit::ModelProjection::Ptr projection = model->makeProjection(psf, wcs, footprint);
    BOOST_CHECK_EQUAL(projection->getModel(), model);
    BOOST_CHECK_EQUAL(projection->getPsfParameterSize(), 0);
    BOOST_CHECK_EQUAL(projection->getWcsParameterSize(), 0);

    multifit::FourierModelProjection::Ptr asFourier = boost::static_pointer_cast<multifit::FourierModelProjection>(projection);
    BOOST_CHECK(asFourier);
    BOOST_CHECK_NO_THROW(model->setLinearParameters(linear.data()));
    BOOST_CHECK_NO_THROW(model->setNonlinearParameters(nonlinear.data()));
}
