// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
#include <array>

#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/geom.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/ellipses/GridTransform.h"
#include "lsst/meas/modelfit/DoubleShapeletPsfApprox.h"

namespace lsst { namespace meas { namespace modelfit {
namespace {
base::FlagDefinitionList flagDefinitions;
} // end anonymous

base::FlagDefinition const DoubleShapeletPsfApproxAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
base::FlagDefinition const DoubleShapeletPsfApproxAlgorithm::INVALID_POINT_FOR_PSF = flagDefinitions.add("flag_invalidPointForPsf", "PSF model could not be evaluated at the source position");
base::FlagDefinition const DoubleShapeletPsfApproxAlgorithm::INVALID_MOMENTS = flagDefinitions.add("flag_invalidMoments", "Moments of the PSF model were not a valid ellipse");
base::FlagDefinition const DoubleShapeletPsfApproxAlgorithm::MAX_ITERATIONS = flagDefinitions.add("flag_maxIterations", "optimizer exceeded the maximum number (inner or outer) iterations");

base::FlagDefinitionList const & DoubleShapeletPsfApproxAlgorithm::getFlagDefinitions() {
    return flagDefinitions;
}


namespace {

// A function that applies a functor to each pixel in an image, passing
// the pixel value and PARENT position to each call.
template <typename Image, typename Func>
void applyPixelFunctor(Image const & image, Func & f) {
    int const height = image.getHeight();
    for (int row = 0; row < height; ++row) {
        int const y = row + image.getY0();
        int x = image.getX0();
        auto const end = image.row_end(row);
        for (auto iter = image.row_begin(row); iter != end; ++iter, ++x) {
            f(*iter, x, y);
        }
    }
}

// A functor for use with applyPixelFunctor that computes unweighted moments.
class MomentsFunctor {
public:

    MomentsFunctor() : _m0(0.0), _mx(0.0), _my(0.0), _mxx(0.0), _myy(0.0), _mxy(0.0) {}

    void operator()(Scalar value, Scalar x, Scalar y) {
        _m0 += value;
        _mx += x*value;
        _my += y*value;
        _mxx += x*x*value;
        _myy += y*y*value;
        _mxy += x*y*value;
    }

    Scalar getSum() const { return _m0; }

    afw::geom::ellipses::Ellipse getEllipse() const {
        geom::Point2D center(_mx/_m0, _my/_m0);
        afw::geom::ellipses::Quadrupole quadrupole(
            _mxx/_m0 - center.getX()*center.getX(),
            _myy/_m0 - center.getY()*center.getY(),
            _mxy/_m0 - center.getX()*center.getY()
        );
        return afw::geom::ellipses::Ellipse(quadrupole, center);
    }

private:
    Scalar _m0;
    Scalar _mx;
    Scalar _my;
    Scalar _mxx;
    Scalar _myy;
    Scalar _mxy;
};

// A functor for use with applyPixelFunctor that flattens the coordinates and
// pixel values into 1-d arrays.
class FlattenFunctor {
public:

    FlattenFunctor(
        ndarray::Array<Scalar,1,1> const & x,
        ndarray::Array<Scalar,1,1> const & y,
        ndarray::Array<Scalar,1,1> const & data
    ) : _xIter(x.begin()), _yIter(y.begin()), _dataIter(data.begin()) {}

    void operator()(Scalar value, Scalar x, Scalar y) {
        *_xIter = x;
        *_yIter = y;
        *_dataIter = value;
        ++_xIter;
        ++_yIter;
        ++_dataIter;
    }

private:
    ndarray::Array<Scalar,1,1>::Iterator _xIter;
    ndarray::Array<Scalar,1,1>::Iterator _yIter;
    ndarray::Array<Scalar,1,1>::Iterator _dataIter;
};

// A functor for use with applyPixelFunctor that flattens data and computes
// the exponential argument (z in A*exp(z)) of a Gaussian function at each point.
class GaussianArgFunctor {
public:

    GaussianArgFunctor(
        ndarray::Array<Scalar,1,1> const & data,
        ndarray::Array<Scalar,1,1> const & arg,
        geom::AffineTransform const * gt
    ) : _dataIter(data.begin()), _argIter(arg.begin()), _gt(gt) {}

    void operator()(Scalar value, Scalar x, Scalar y) {
        *_dataIter = value;
        *_argIter = -0.5*(*_gt)(geom::Point2D(x, y)).asEigen().squaredNorm();
        ++_dataIter;
        ++_argIter;
    }

private:
    ndarray::Array<Scalar,1,1>::Iterator _dataIter;
    ndarray::Array<Scalar,1,1>::Iterator _argIter;
    geom::AffineTransform const * _gt;
};

} // anonymous


DoubleShapeletPsfApproxAlgorithm::DoubleShapeletPsfApproxAlgorithm(
    DoubleShapeletPsfApproxControl const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _centroidExtractor(schema, name)
{
    std::vector<int> const orders = { ctrl.innerOrder, ctrl.outerOrder };
    _key = shapelet::MultiShapeletFunctionKey::addFields(
        schema,
        name,
        "Double-Shapelet approximation to the PSF model at the position of this source",
        "pixel",  // ellipse units
        "",       // amplitude is unitless
        orders
    );
    _flagHandler = meas::base::FlagHandler::addFields(schema, name, getFlagDefinitions());
}


shapelet::MultiShapeletFunction DoubleShapeletPsfApproxAlgorithm::initializeResult(Control const & ctrl) {
    // Construct inner component with unit flux, unit circle ellipse.
    shapelet::ShapeletFunction inner(ctrl.innerOrder, shapelet::HERMITE);
    inner.getCoefficients()[0] = 1.0 / shapelet::ShapeletFunction::FLUX_FACTOR;
    // Construct outer component with unit flux, circle with radius = ctrl.radiusRatio.
    shapelet::ShapeletFunction outer(ctrl.outerOrder, shapelet::HERMITE, ctrl.radiusRatio);
    outer.getCoefficients()[0] = 1.0 / shapelet::ShapeletFunction::FLUX_FACTOR;
    // Rescale outer component's amplitude to get peak ratio right.
    double currentPeakRatio =  outer.evaluate()(0.0, 0.0) / inner.evaluate()(0.0, 0.0);
    outer.getCoefficients()[0] *= ctrl.peakRatio / currentPeakRatio;
    // Construct the compound object.
    shapelet::MultiShapeletFunction result;
    result.getComponents().push_back(std::move(inner));
    result.getComponents().push_back(std::move(outer));
    // Normalize to unit flux and unit circle moments.
    result.normalize();
    result.transformInPlace(result.evaluate().computeMoments().getGridTransform());
    return result;
}

void DoubleShapeletPsfApproxAlgorithm::fitMoments(
    shapelet::MultiShapeletFunction & result,
    Control const & ctrl,
    afw::detection::Psf::Image const & psfImage
) {
    MomentsFunctor func;
    applyPixelFunctor(psfImage, func);
    auto ellipse = func.getEllipse();
    auto sum = func.getSum();
    try {
        ellipse.getCore().normalize();
    } catch (pex::exceptions::InvalidParameterError & err) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            err.what(),
            INVALID_MOMENTS.number
        );
    }
    result.transformInPlace(ellipse.getGridTransform().inverted());
    result.normalize(sum);
    afw::geom::ellipses::Axes axesInner(result.getComponents()[0].getEllipse().getCore());
    afw::geom::ellipses::Axes axesOuter(result.getComponents()[1].getEllipse().getCore());
    if (!(axesInner.getB() >= ctrl.minRadius)) {  // compare for opposite to catch NaNs
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            (boost::format("Semi-minor radius derived from moments (%g) is less than lower bound (%g).")
             % axesOuter.getB() % ctrl.minRadius).str(),
            INVALID_MOMENTS.number
        );
    }
    double const maxRadius = ctrl.maxRadiusBoxFraction*std::sqrt(psfImage.getBBox().getArea());
    if (!(axesOuter.getA() <= maxRadius)) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            (boost::format("Semi-major radius derived from moments (%g) is greater than upper bound (%g).")
             % axesOuter.getA() % maxRadius).str(),
            INVALID_MOMENTS.number
        );
    }
    if (!(axesOuter.getDeterminantRadius() - axesInner.getDeterminantRadius() >= ctrl.minRadiusDiff)) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            (boost::format("Different between radii (%g - %g = %g) is smaller than lower bound (%g).")
             % axesOuter.getDeterminantRadius() % axesInner.getDeterminantRadius()
             % (axesOuter.getDeterminantRadius() - axesInner.getDeterminantRadius())
             % ctrl.minRadiusDiff).str(),
            INVALID_MOMENTS.number
        );
    }
}

namespace {

// An OptimizerObjective that fits the profile of a DoubleShapelet model
// (see DoubleShapeletPsfApproxAlgorithm::makeObjective for more info.)
class ProfileObjective : public OptimizerObjective {
public:

    ProfileObjective(
        afw::geom::ellipses::Ellipse const & moments,
        DoubleShapeletPsfApproxControl const & ctrl,
        afw::image::Image<Scalar> const & image
    ) : OptimizerObjective(image.getBBox().getArea(), 4),
        _minRadius(ctrl.minRadius),
        _maxRadius(ctrl.maxRadiusBoxFraction * std::sqrt(this->dataSize)),
        _minRadiusDiff(ctrl.minRadiusDiff),
        _data(ndarray::Array<Scalar,1,1>(ndarray::allocate(this->dataSize))),
        _arg(ndarray::Array<Scalar,1,1>(ndarray::allocate(this->dataSize)))
    {
        // Radius parameters are defined as factors of the moments ellipse, so
        // we have to scale the constraints the same way.
        afw::geom::ellipses::Axes axes(moments.getCore());
        _minRadius /= axes.getB();
        _maxRadius /= axes.getA();
        _minRadiusDiff /= axes.getDeterminantRadius();
        // We now compute most of the exponential argument to the Gaussian function up
        // front (everything but the factor of 1/r^2), since that won't change as
        // the parameters change.
        geom::AffineTransform gt = moments.getGridTransform();
        GaussianArgFunctor func(_data.shallow(), _arg.shallow(), &gt);
        applyPixelFunctor(image, func);
        _normalization = gt.getLinear().computeDeterminant() * shapelet::ShapeletFunction::FLUX_FACTOR /
            (2.0*geom::PI);
    }

    virtual void computeResiduals(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & residuals
    ) const {
        Scalar iAlpha = parameters[0] * _normalization;
        Scalar oAlpha = parameters[1] * _normalization;
        Scalar iR2 = parameters[2] * parameters[2];
        Scalar oR2 = parameters[3] * parameters[3];
        auto argEigen = ndarray::asEigenArray(_arg);
        ndarray::asEigenArray(residuals) = ndarray::asEigenArray(_data)
            - (iAlpha/iR2)*(argEigen/iR2).exp()
            - (oAlpha/oR2)*(argEigen/oR2).exp();
    }

    virtual bool differentiateResiduals(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,2,-2> const & derivatives
    ) const {
        auto d = ndarray::asEigenArray(derivatives);
        auto argEigen = ndarray::asEigenArray(_arg);
        Scalar iR2 = parameters[2] * parameters[2];
        Scalar oR2 = parameters[3] * parameters[3];
        d.col(0) = - (_normalization/iR2)*(argEigen/iR2).exp();
        d.col(1) = - (_normalization/oR2)*(argEigen/oR2).exp();
        d.col(2) = -2.0*d.col(0)*(argEigen/iR2 + 1.0)*parameters[0]/parameters[2];
        d.col(3) = -2.0*d.col(1)*(argEigen/oR2 + 1.0)*parameters[1]/parameters[3];
        return true;
    }

    virtual bool hasPrior() const { return true; }

    virtual Scalar computePrior(ndarray::Array<Scalar const,1,1> const & parameters) const {
        return ((parameters[3] - parameters[2] < _minRadiusDiff)
                || (parameters[3] > _maxRadius)
                || (parameters[2] < _minRadius))
            ? 0.0 : 1.0;
    }

    virtual void differentiatePrior(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,1> const & hessian
    ) const {
        gradient.deep() = 0.0;
        hessian.deep() = 0.0;
    }

private:
    Scalar _minRadius;
    Scalar _maxRadius;
    Scalar _minRadiusDiff;
    Scalar _normalization;
    ndarray::Array<Scalar,1,1> _data;
    ndarray::Array<Scalar,1,1> _arg;
};

} // anonymous

std::shared_ptr<OptimizerObjective> DoubleShapeletPsfApproxAlgorithm::makeObjective(
    afw::geom::ellipses::Ellipse const & moments,
    Control const & ctrl,
    afw::image::Image<Scalar> const & psfImage
) {
    return std::make_shared<ProfileObjective>(moments, ctrl, psfImage);
}

void DoubleShapeletPsfApproxAlgorithm::fitProfile(
    shapelet::MultiShapeletFunction & result,
    Control const & ctrl,
    afw::detection::Psf::Image const & psfImage
) {
    afw::geom::ellipses::Ellipse moments = result.evaluate().computeMoments();
    Scalar momentsRadius = moments.getCore().getDeterminantRadius();
    auto objective = makeObjective(moments, ctrl, psfImage);
    ndarray::Array<Scalar,1,1> parameters = ndarray::allocate(objective->parameterSize);
    parameters[0] = result.getComponents()[0].getCoefficients()[0];
    parameters[1] = result.getComponents()[1].getCoefficients()[0];
    parameters[2] = result.getComponents()[0].getEllipse().getCore().getDeterminantRadius() / momentsRadius;
    parameters[3] = result.getComponents()[1].getEllipse().getCore().getDeterminantRadius() / momentsRadius;
    Optimizer optimizer(objective, parameters, ctrl.optimizer);
    optimizer.run();
    int state = optimizer.getState();
    result.getComponents()[0].getCoefficients()[0] = optimizer.getParameters()[0];
    result.getComponents()[1].getCoefficients()[0] = optimizer.getParameters()[1];
    result.getComponents()[0].getEllipse().getCore().scale(optimizer.getParameters()[2] / parameters[2]);
    result.getComponents()[1].getEllipse().getCore().scale(optimizer.getParameters()[3] / parameters[3]);
    if (state & Optimizer::FAILED) {
        if (state & Optimizer::FAILED_MAX_ITERATIONS) {
            throw LSST_EXCEPT(
                meas::base::MeasurementError,
                "Fitting profile failed; too many outer iterations",
                MAX_ITERATIONS.number
            );
        } else {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                (boost::format(
                    "Unexpected error (state=%#04x) while fitting profile; please report this as a bug"
                 ) % state).str()
            );
        }
    }
}

void DoubleShapeletPsfApproxAlgorithm::fitShapelets(
    shapelet::MultiShapeletFunction & result,
    Control const & ctrl,
    afw::detection::Psf::Image const & psfImage
) {
    if (result.getComponents().size() != 2u) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "MultiShapeletFunction argment to fitShapelets must have exactly 2 components."
        );
    }
    if (ctrl.innerOrder == 0 && ctrl.outerOrder == 0) {
        // No higher-order shapelet terms to fit.
        return;
    }
    auto & innerComponent = result.getComponents().front();
    auto & outerComponent = result.getComponents().back();
    // Create flattened coordinate arrays to pass to MatrixBuilders, while copying pixel values into
    // another flattened array.
    int area = psfImage.getBBox().getArea();
    ndarray::Array<Scalar,1,1> xArray = ndarray::allocate(area);
    ndarray::Array<Scalar,1,1> yArray = ndarray::allocate(area);
    ndarray::Array<Scalar,1,1> data = ndarray::allocate(area);
    FlattenFunctor func(xArray, yArray, data);
    applyPixelFunctor(psfImage, func);
    // Construct two MatrixBuilders, using the pattern that lets them share workspace.
    shapelet::MatrixBuilder<Scalar> innerBuilder(xArray, yArray, ctrl.innerOrder);
    ndarray::Array<Scalar,2,-2> innerMatrix = innerBuilder(innerComponent.getEllipse());
    int n1 = innerBuilder.getBasisSize();
    shapelet::MatrixBuilder<Scalar> outerBuilder(xArray, yArray, ctrl.outerOrder);
    ndarray::Array<Scalar,2,-2> outerMatrix = outerBuilder(outerComponent.getEllipse());
    int n2 = outerBuilder.getBasisSize();
    // Build the full matrix we'll need to solve: the combination of the two basis, without
    // zeroth-order terms.
    int fitBasisSize = n1 + n2 - 2;
    ndarray::Array<Scalar,2,-2> fitMatrix = ndarray::allocate(data.size(), fitBasisSize);
    if (n1 > 1) {
        fitMatrix[ndarray::view()(0, n1 - 1)] = innerMatrix[ndarray::view()(1, n1)];
    }
    if (n2 > 1) {
        fitMatrix[ndarray::view()(n1 - 1, n1 + n2 - 2)] = outerMatrix[ndarray::view()(1, n2)];
    }
    // Subtract the zeroth-order model from the data.
    ndarray::asEigenMatrix(data) -=
            ndarray::asEigenMatrix(innerMatrix).col(0) * innerComponent.getCoefficients()[0];
    ndarray::asEigenMatrix(data) -=
            ndarray::asEigenMatrix(outerMatrix).col(0) * outerComponent.getCoefficients()[0];
    // Finaly, we can do the fit, and stuff the fitted coefficients back into the result.
    auto lstsq = afw::math::LeastSquares::fromDesignMatrix(fitMatrix, data);
    if (n1 > 1) {
        innerComponent.getCoefficients()[ndarray::view(1, n1)] =
            lstsq.getSolution()[ndarray::view(0, n1 - 1)];
    }
    if (n2 > 1) {
        outerComponent.getCoefficients()[ndarray::view(1, n2)] =
            lstsq.getSolution()[ndarray::view(n1 - 1, n1 + n2 - 2)];
    }
}

void DoubleShapeletPsfApproxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    auto psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            meas::base::FatalAlgorithmError,
            "No Psf attached to Exposure for DoubleShapeletPsfApprox."
        );
    }
    auto position = _centroidExtractor(measRecord, _flagHandler);
    std::shared_ptr<afw::detection::Psf::Image> psfImage;
    try {
        psfImage = psf->computeKernelImage(position);
    } catch (pex::exceptions::Exception & err) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            err.what(),
            INVALID_POINT_FOR_PSF.number
        );
    }
    auto result = initializeResult(_ctrl);
    fitMoments(result, _ctrl, *psfImage);
    fitProfile(result, _ctrl, *psfImage);
    fitShapelets(result, _ctrl, *psfImage);
    measRecord.set(_key, result);
}


void DoubleShapeletPsfApproxAlgorithm::fail(
    afw::table::SourceRecord & measRecord,
    lsst::meas::base::MeasurementError * error
) const {
    _flagHandler.handleFailure(measRecord, error);
}


}}} // namespace lsst::meas::modelfit
