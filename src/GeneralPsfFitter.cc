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

#include "ndarray/eigen.h"

#include "lsst/pex/exceptions.h"
#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/modelfit/GeneralPsfFitter.h"

namespace lsst { namespace meas { namespace modelfit {
namespace {
base::FlagDefinitionList flagDefinitions;
} // end anonymous

base::FlagDefinition const GeneralPsfFitterAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
base::FlagDefinition const GeneralPsfFitterAlgorithm::MAX_INNER_ITERATIONS = flagDefinitions.add("flag_max_inner_iterations", "exceeded maxInnerIterations");
base::FlagDefinition const GeneralPsfFitterAlgorithm::MAX_OUTER_ITERATIONS = flagDefinitions.add("flag_max_outer_iterations", "exceeded maxOuterIterations");
base::FlagDefinition const GeneralPsfFitterAlgorithm::EXCEPTION = flagDefinitions.add("flag_exception", "exception in apply method");
base::FlagDefinition const GeneralPsfFitterAlgorithm::CONTAINS_NAN = flagDefinitions.add("flag_contains_nan", "Nan in the Psf image");

base::FlagDefinitionList const & GeneralPsfFitterAlgorithm::getFlagDefinitions() {
    return flagDefinitions;
}


namespace {

typedef std::pair<std::string,GeneralPsfFitterComponentControl> Component;
typedef std::vector<Component> ComponentVector;
typedef ComponentVector::const_iterator ComponentIterator;

struct ComponentNameIs {

    explicit ComponentNameIs(std::string const & target) : _target(target) {}

    bool operator()(Component const & a) const {
        return a.first == _target;
    }

private:
    std::string _target;
};

static afw::geom::ellipses::SeparableConformalShearLogTraceRadius psfFitterEllipseCore;

class GeneralPsfFitterModel : public Model {
public:

    GeneralPsfFitterModel(
        BasisVector basisVector,
        NameVector nonlinearNames,
        NameVector amplitudeNames,
        NameVector fixedNames,
        ComponentVector components
    ) : Model(basisVector, nonlinearNames, amplitudeNames, fixedNames), _components(components) {}

    virtual std::shared_ptr<Prior> adaptPrior(std::shared_ptr<Prior> prior) const {
        if (prior->getTag() != "PSF") {
            throw LSST_EXCEPT(pex::exceptions::LogicError, "Invalid prior for this model");
        }
        return prior;
    }

    virtual EllipseVector makeEllipseVector() const {
        return EllipseVector(getBasisCount(), afw::geom::ellipses::Ellipse(psfFitterEllipseCore));
    }

    virtual void writeEllipses(
        Scalar const * nonlinearIter, Scalar const * fixedIter, EllipseIterator ellipseIter
    ) const {
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i, ++ellipseIter) {
            afw::geom::ellipses::Ellipse::ParameterVector v;
            // we always store fiducial values in the fixed parameter array; if a parameter is allowed
            // to vary, we store the offset from the fiducial value in the nonlinear parameter array.
            for (int n = 0; n < 5; ++n, ++fixedIter) {
                v[n] = *fixedIter;
            }
            if (i->second.ellipticityPriorSigma > 0.0) {
                v[0] += *(nonlinearIter++);
                v[1] += *(nonlinearIter++);
            }
            if (i->second.radiusPriorSigma > 0.0) {
                v[2] += *(nonlinearIter++);
            }
            if (i->second.positionPriorSigma > 0.0) {
                v[3] += *(nonlinearIter++);
                v[4] += *(nonlinearIter++);
            }
            ellipseIter->setParameterVector(v);
        }
    }

    virtual void readEllipses(
        EllipseConstIterator ellipseIter, Scalar * nonlinearIter, Scalar * fixedIter
    ) const {
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i, ++ellipseIter) {
            afw::geom::ellipses::Ellipse::ParameterVector v = ellipseIter->getParameterVector();
            // We always store fiducial values in the fixed parameter array; if a parameter is allowed
            // to vary, we store the offset from the fiducial value in the nonlinear parameter array.
            // When reading ellipses, we only update the fixed array, and set the offsets in the
            // nonlinear array to zero.
            // I'm a little concerned that this means we can't round-trip parameter vectors through
            // ellipses with this model, but I can't think of a concrete case where that would cause
            // problems.
            for (int n = 0; n < 5; ++n, ++fixedIter) {
                *fixedIter = v[n];
            }
            if (i->second.ellipticityPriorSigma > 0.0) {
                *(nonlinearIter++) = 0;
                *(nonlinearIter++) = 0;
            }
            if (i->second.radiusPriorSigma > 0.0) {
                *(nonlinearIter++) = 0;
            }
            if (i->second.positionPriorSigma > 0.0) {
                *(nonlinearIter++) = 0;
                *(nonlinearIter++) = 0;
            }
        }
    }

    std::size_t findComponentIndex(std::string const & name) const {
        ComponentIterator i = std::find_if(_components.begin(), _components.end(), ComponentNameIs(name));
        return i - _components.begin();
    }

    shapelet::MultiShapeletFunction adapt(
        shapelet::MultiShapeletFunction const & previousFit,
        GeneralPsfFitterModel const & previousModel
    ) const {
        shapelet::MultiShapeletFunction result;
        shapelet::ShapeletFunction const & previousPrimary
            = previousFit.getComponents()[previousModel.findComponentIndex("primary")];
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i) {
            result.getComponents().push_back(shapelet::ShapeletFunction(i->second.order, shapelet::HERMITE));
            shapelet::ShapeletFunction & current = result.getComponents().back();
            std::size_t n = previousModel.findComponentIndex(i->first);
            if (n >= previousFit.getComponents().size()) {
                // previous didn't include the component we're setting up now, so we base its parameters
                // off the previous primary component (but we initialize the coefficients to zero).
                current.setEllipse(previousPrimary.getEllipse());
                current.getEllipse().getCore().scale(i->second.radiusFactor);
            } else {
                // previous did include the component we're setting up now, so we copy it over,
                // including as many coeffients as possible.
                shapelet::ShapeletFunction const & previousComponent = previousFit.getComponents()[n];
                current.setEllipse(previousComponent.getEllipse());
                int minOrder = std::min(current.getOrder(), previousComponent.getOrder());
                int minSize = shapelet::computeSize(minOrder);
                current.getCoefficients()[ndarray::view(0, minSize)]
                    = previousComponent.getCoefficients()[ndarray::view(0, minSize)];
            }
        }
        return result;
    }

    shapelet::MultiShapeletFunction makeInitial(afw::geom::ellipses::Quadrupole const & moments) const {
        shapelet::MultiShapeletFunction result;
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i) {
            result.getComponents().push_back(shapelet::ShapeletFunction(i->second.order, shapelet::HERMITE));
            shapelet::ShapeletFunction & current = result.getComponents().back();
            current.getEllipse().setCore(moments);
            current.getEllipse().getCore().scale(i->second.radiusFactor);
            if (i->first == "primary") {
                current.getCoefficients()[0] = 1.0 / shapelet::ShapeletFunction::FLUX_FACTOR;
            }
        }
        return result;
    }

    void fillParameters(
        shapelet::MultiShapeletFunction const & input,
        ndarray::Array<Scalar,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & fixed
    ) const {
        EllipseVector ellipses = makeEllipseVector();
        LSST_THROW_IF_NE(
            input.getComponents().size(), ellipses.size(),
            pex::exceptions::LengthError,
            "Number of elements in input multishapelet (%d) does not match expected number of elements (%d)"
        );
        assert(ellipses.size() == _components.size()); // should be guaranteed by construction
        ndarray::Array<Scalar,1,1>::Iterator amplitudeIter = amplitudes.begin();
        for (std::size_t n = 0; n < ellipses.size(); ++n) {
            shapelet::ShapeletFunction const & current = input.getComponents()[n];
            if (current.getOrder() != _components[n].second.order) {
                throw LSST_EXCEPT(
                    pex::exceptions::LengthError,
                    (boost::format("Shapelet order of component %d has order %d, not %d")
                     % n % current.getOrder() % _components[n].second.order).str()
                );
            }
            ellipses[n] = current.getEllipse();
            std::copy(current.getCoefficients().begin(), current.getCoefficients().end(), amplitudeIter);
            amplitudeIter += current.getCoefficients().getSize<0>();
        }
        readEllipses(ellipses.begin(), nonlinear.begin(), fixed.begin());
    }

private:
    ComponentVector _components;
};

class GeneralPsfFitterPrior : public Prior {
public:

    GeneralPsfFitterPrior(ComponentVector components) : Prior("PSF"), _components(components) {}

    Scalar square(Scalar x) const { return x*x; }

    virtual Scalar evaluate(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes
    ) const {
        Scalar z = 0;
        ndarray::Array<Scalar const,1,1>::Iterator nonlinearIter = nonlinear.begin();
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i) {
            // If we have a non-delta-fn prior, it's gaussian with zero mean in the nonlinear parameter,
            // as we've defined that as the offset between the actual ellipse parameter and the
            // fiducial value.
            if (i->second.ellipticityPriorSigma > 0.0) {
                z += square(*(nonlinearIter++) / i->second.ellipticityPriorSigma);
                z += square(*(nonlinearIter++) / i->second.ellipticityPriorSigma);
            }
            if (i->second.radiusPriorSigma > 0.0) {
                z += square(*(nonlinearIter++) / i->second.radiusPriorSigma);
            }
            if (i->second.positionPriorSigma > 0.0) {
                z += square(*(nonlinearIter++) / i->second.positionPriorSigma);
                z += square(*(nonlinearIter++) / i->second.positionPriorSigma);
            }
        }
        assert(nonlinearIter == nonlinear.end());
        return std::exp(-0.5*z);
    }

    virtual void evaluateDerivatives(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & nonlinearGradient,
        ndarray::Array<Scalar,1,1> const & amplitudeGradient,
        ndarray::Array<Scalar,2,1> const & nonlinearHessian,
        ndarray::Array<Scalar,2,1> const & amplitudeHessian,
        ndarray::Array<Scalar,2,1> const & crossHessian
    ) const {
        Scalar p = evaluate(nonlinear, amplitudes);
        int n = 0;
        nonlinearHessian.deep() = 0.0;
        for (ComponentIterator i = _components.begin(); i != _components.end(); ++i) {
            if (i->second.ellipticityPriorSigma > 0.0) {
                Scalar sigmaSqrInv = square(1.0/i->second.ellipticityPriorSigma);
                for (int k = 0; k < 2; ++k, ++n) {
                    nonlinearGradient[n] = -nonlinear[n]*sigmaSqrInv;
                    nonlinearHessian[n][n] = -sigmaSqrInv;
                }
            }
            if (i->second.radiusPriorSigma > 0.0) {
                Scalar sigmaSqrInv = square(1.0/i->second.radiusPriorSigma);
                for (int k = 0; k < 1; ++k, ++n) {
                    nonlinearGradient[n] = -nonlinear[n]*sigmaSqrInv;
                    nonlinearHessian[n][n] = -sigmaSqrInv;
                }
            }
            if (i->second.positionPriorSigma > 0.0) {
                Scalar sigmaSqrInv = square(1.0/i->second.positionPriorSigma);
                for (int k = 0; k < 2; ++k, ++n) {
                    nonlinearGradient[n] = -nonlinear[n]*sigmaSqrInv;
                    nonlinearHessian[n][n] = -sigmaSqrInv;
                }
            }
        }
        auto nonlinearHessianMatrixMap = ndarray::asEigenMatrix(nonlinearHessian);
        auto nonlinearGradientMatrixMap = ndarray::asEigenMatrix(nonlinearGradient);
        nonlinearHessianMatrixMap.selfadjointView<Eigen::Lower>().rankUpdate(nonlinearGradientMatrixMap);
        nonlinearHessianMatrixMap = nonlinearHessianMatrixMap.selfadjointView<Eigen::Lower>();
        nonlinearGradientMatrixMap *= p;
        nonlinearHessianMatrixMap *= p;
        amplitudeGradient.deep() = 0.0;
        amplitudeHessian.deep() = 0.0;
        crossHessian.deep() = 0.0;
    }

    virtual Scalar marginalize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear
    ) const {
        // Don't need this unless we want to sample PSF models
        throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
    }

    virtual Scalar maximize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & amplitudes
    ) const {
        // Don't need this unless we want to sample PSF models
        throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
    }

    virtual void drawAmplitudes(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        afw::math::Random & rng,
        ndarray::Array<Scalar,2,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & weights,
        bool multiplyWeights=false
    ) const {
        // Don't need this unless we want to sample PSF models
        throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
    }

private:
    ComponentVector _components;
};

ComponentVector vectorizeComponents(GeneralPsfFitterControl const & ctrl) {
    ComponentVector components;
    if (ctrl.inner.order >= 0) {
        components.push_back(Component("inner", ctrl.inner));
    }
    if (ctrl.primary.order >= 0) {
        components.push_back(Component("primary", ctrl.primary));
    }
    if (ctrl.wings.order >= 0) {
        components.push_back(Component("wings", ctrl.wings));
    }
    if (ctrl.outer.order >= 0) {
        components.push_back(Component("outer", ctrl.outer));
    }
    return components;
}

} // anonymous

GeneralPsfFitter::GeneralPsfFitter(GeneralPsfFitterControl const & ctrl) :
    _ctrl(ctrl)
{
    if (_ctrl.primary.order < 0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "GeneralPsfFitter control must have a primary component with nonnegative order"
        );
    }
    ComponentVector components = vectorizeComponents(_ctrl);

    _prior = std::make_shared<GeneralPsfFitterPrior>(components);

    Model::BasisVector basisVector;
    Model::NameVector nonlinearNames;
    Model::NameVector amplitudeNames;
    Model::NameVector fixedNames;

    for (ComponentIterator i = components.begin(); i != components.end(); ++i) {

        // Construct a MultiShapeletBasis with a single shapelet basis with this component
        int dim = shapelet::computeSize(i->second.order);
        std::shared_ptr<shapelet::MultiShapeletBasis> basis = std::make_shared<shapelet::MultiShapeletBasis>(dim);
        ndarray::Array<double,2,2> matrix = ndarray::allocate(dim, dim);
        ndarray::asEigenMatrix(matrix).setIdentity();
        basis->addComponent(1.0, i->second.order, matrix);
        basisVector.push_back(basis);

        // Append to the name vectors for this component; all of this has to be consistent with the
        // iteration that happens in GeneralPsfFitterModel and GeneralPsfFitterPrior
        for (shapelet::PackedIndex s; s.getOrder() <= i->second.order; ++s) {
            amplitudeNames.push_back(
                (boost::format("%s.alpha[%d,%d]") % i->first % s.getX() % s.getY()).str()
            );
        }
        fixedNames.push_back(i->first + ".fiducial.eta1");
        fixedNames.push_back(i->first + ".fiducial.eta2");
        fixedNames.push_back(i->first + ".fiducial.logR");
        fixedNames.push_back(i->first + ".fiducial.x");
        fixedNames.push_back(i->first + ".fiducial.y");
        if (i->second.ellipticityPriorSigma > 0.0) {
            nonlinearNames.push_back(i->first + ".eta1");
            nonlinearNames.push_back(i->first + ".eta2");
        }
        if (i->second.radiusPriorSigma > 0.0) {
            nonlinearNames.push_back(i->first + ".logR");
        }
        if (i->second.positionPriorSigma > 0.0) {
            nonlinearNames.push_back(i->first + ".x");
            nonlinearNames.push_back(i->first + ".y");
        }

    }

    _model = std::make_shared<GeneralPsfFitterModel>(
        basisVector, nonlinearNames, amplitudeNames, fixedNames, components
    );
}

shapelet::MultiShapeletFunctionKey GeneralPsfFitter::addFields(
    afw::table::Schema & schema,
    std::string const & prefix
) const {
    ComponentVector components = vectorizeComponents(_ctrl);
    std::vector<int> orders(components.size());
    for (std::size_t i = 0; i < components.size(); ++i) {
        orders[i] = components[i].second.order;
    }
    return shapelet::MultiShapeletFunctionKey::addFields(
        schema, prefix, "multi-Shapelet approximation to the PSF model",
        "pixel", // ellipse units
        "",       // coefficient units (unitless)
        orders
    );
}

shapelet::MultiShapeletFunction GeneralPsfFitter::adapt(
    shapelet::MultiShapeletFunction const & previousFit,
    std::shared_ptr<Model> previousModel
) const {
    std::shared_ptr<GeneralPsfFitterModel> m = std::dynamic_pointer_cast<GeneralPsfFitterModel>(previousModel);
    if (!m) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "Model passed to GeneralPsfFitter::adapt must have been constructed by GeneralPsfFitter"
        );
    }
    return std::static_pointer_cast<GeneralPsfFitterModel>(_model)->adapt(previousFit, *m);
}


shapelet::MultiShapeletFunction GeneralPsfFitter::apply(
    afw::image::Image<Pixel> const & image,
    afw::geom::ellipses::Quadrupole const & moments,
    Scalar noiseSigma,
    int * pState
) const {
    if (noiseSigma <= 0) {
        noiseSigma = _ctrl.defaultNoiseSigma;
    }
    shapelet::MultiShapeletFunction initial
        = std::static_pointer_cast<GeneralPsfFitterModel>(_model)->makeInitial(moments);
    return apply(image, initial, noiseSigma, pState);
}

shapelet::MultiShapeletFunction GeneralPsfFitter::apply(
    afw::image::Image<Pixel> const & image,
    shapelet::MultiShapeletFunction const & initial,
    Scalar noiseSigma,
    int * pState
) const {
    if (noiseSigma <= 0) {
        noiseSigma = _ctrl.defaultNoiseSigma;
    }
    int const parameterDim = _model->getNonlinearDim() + _model->getAmplitudeDim();
    ndarray::Array<Scalar,1,1> parameters = ndarray::allocate(parameterDim);
    ndarray::Array<Scalar,1,1> nonlinear = parameters[ndarray::view(0, _model->getNonlinearDim())];
    ndarray::Array<Scalar,1,1> amplitudes
        = parameters[ndarray::view(_model->getNonlinearDim(), parameterDim)];
    ndarray::Array<Scalar,1,1> fixed = ndarray::allocate(_model->getFixedDim());

    std::static_pointer_cast<GeneralPsfFitterModel>(_model)->fillParameters(initial, nonlinear, amplitudes, fixed);

    std::shared_ptr<Likelihood> likelihood = std::make_shared<MultiShapeletPsfLikelihood>(
        image.getArray(), image.getXY0(), _model, noiseSigma, fixed
    );
    std::shared_ptr<OptimizerObjective> objective = OptimizerObjective::makeFromLikelihood(likelihood, _prior);
    Optimizer optimizer(objective, parameters, _ctrl.optimizer);
    optimizer.run();

    parameters.deep() = optimizer.getParameters(); // this sets nonlinear, amplitudes, because they're views
    if (pState != nullptr) {
        *pState = optimizer.getState();
    }
    return _model->makeShapeletFunction(nonlinear, amplitudes, fixed);
}

GeneralPsfFitterAlgorithm::GeneralPsfFitterAlgorithm(GeneralPsfFitterControl const & ctrl,
    afw::table::Schema & schema,
    std::string const & prefix
) : GeneralPsfFitter(ctrl)
{
    _flagHandler = lsst::meas::base::FlagHandler::addFields(schema, prefix,
                                          getFlagDefinitions());
    _key = addFields(schema, prefix);
}

void GeneralPsfFitterAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Image<double> const & image,
    shapelet::MultiShapeletFunction const & initial
) const {
    int state = 0;
    shapelet::MultiShapeletFunction result = apply(image, initial, -1, &state);
    measRecord.set(_key, result);
    if (state & Optimizer::FAILED_MAX_INNER_ITERATIONS) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            MAX_INNER_ITERATIONS.doc,
            MAX_INNER_ITERATIONS.number
        );
    }
    if (state & Optimizer::FAILED_MAX_OUTER_ITERATIONS) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            MAX_OUTER_ITERATIONS.doc,
            MAX_OUTER_ITERATIONS.number
        );
    }
    if (state & Optimizer::FAILED_EXCEPTION) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            EXCEPTION.doc,
            EXCEPTION.number
        );
    }
    if (state & Optimizer::FAILED_NAN) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            CONTAINS_NAN.doc,
            CONTAINS_NAN.number
        );
    }
}

void GeneralPsfFitterAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Image<double> const & image,
    afw::geom::ellipses::Quadrupole const & moments
) const {
    int state = 0;
    shapelet::MultiShapeletFunction result = apply(image, moments, -1, &state);
    measRecord.set(_key, result);
    if (state & Optimizer::FAILED_MAX_INNER_ITERATIONS) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            MAX_INNER_ITERATIONS.doc,
            MAX_INNER_ITERATIONS.number
        );
    }
    if (state & Optimizer::FAILED_MAX_OUTER_ITERATIONS) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            MAX_OUTER_ITERATIONS.doc,
            MAX_OUTER_ITERATIONS.number
        );
    }
    if (state & Optimizer::FAILED_EXCEPTION) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            EXCEPTION.doc,
            EXCEPTION.number
        );
    }
    if (state & Optimizer::FAILED_NAN) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            CONTAINS_NAN.doc,
            CONTAINS_NAN.number
        );
    }
}

void GeneralPsfFitterAlgorithm::fail(
    afw::table::SourceRecord & measRecord,
    lsst::meas::base::MeasurementError * error
) const {
   if (error == nullptr) {
       _flagHandler.handleFailure(measRecord);
   } else {
       _flagHandler.handleFailure(measRecord, error);
   }
}


class MultiShapeletPsfLikelihood::Impl {
public:

    explicit Impl(
        ndarray::Array<Pixel const,1,1> const & x,
        ndarray::Array<Pixel const,1,1> const & y,
        Model::EllipseVector const & ellipses,
        Model::BasisVector const & basisVector,
        Scalar sigma
    ) : _ellipses(ellipses),
        _builders(),
        _sigma(sigma)
    {
        FactoryVector factories;
        factories.reserve(basisVector.size());
        _builders.reserve(basisVector.size());
        int workspaceSize = 0;
        for (Model::BasisVector::const_iterator i = basisVector.begin(); i != basisVector.end(); ++i) {
            factories.push_back(shapelet::MatrixBuilderFactory<Pixel>(x, y, **i));
            workspaceSize = std::max(workspaceSize, factories.back().computeWorkspace());
        }
        shapelet::MatrixBuilderWorkspace<Pixel> workspace(workspaceSize);
        for (FactoryVector::const_iterator i = factories.begin(); i != factories.end(); ++i) {
            shapelet::MatrixBuilderWorkspace<Pixel> wsCopy(workspace); // share workspace between builders
            _builders.push_back((*i)(wsCopy));
        }
    }

    void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & fixed,
        Model const & model
    ) {
        model.writeEllipses(nonlinear.begin(), fixed.begin(), _ellipses.begin());
        modelMatrix.deep() = 0.0;
        Model::BasisVector const & basisVector = model.getBasisVector();
        int amplitudeOffset = 0;
        for (std::size_t i = 0; i < basisVector.size(); ++i) {
            int amplitudeEnd = amplitudeOffset + _builders[i].getBasisSize();
            _builders[i](modelMatrix[ndarray::view()(amplitudeOffset, amplitudeEnd)], _ellipses[i]);
            amplitudeOffset = amplitudeEnd;
        }
        ndarray::asEigenMatrix(modelMatrix) /= _sigma;
    }

private:
    typedef std::vector< shapelet::MatrixBuilder<Pixel> > BuilderVector;
    typedef std::vector< shapelet::MatrixBuilderFactory<Pixel> > FactoryVector;

    Model::EllipseVector _ellipses;
    BuilderVector _builders;
    Scalar _sigma;
};

MultiShapeletPsfLikelihood::MultiShapeletPsfLikelihood(
    ndarray::Array<Pixel const,2,1> const & image,
    geom::Point2I const & xy0,
    std::shared_ptr<Model> model,
    Scalar sigma,
    ndarray::Array<Scalar const,1,1> const & fixed
) :
    Likelihood(model, fixed)
{
    int nx = image.getSize<1>();
    int ny = image.getSize<0>();
    int nTot = nx*ny;
    ndarray::Array<Pixel,1,1> x = ndarray::allocate(nTot);
    ndarray::Array<Pixel,1,1> y = ndarray::allocate(nTot);
    int j = 0;
    for (int iy = xy0.getY(), yEnd = xy0.getY() + ny; iy < yEnd; ++iy) {
        for (int ix = xy0.getX(), xEnd = xy0.getX() + nx; ix < xEnd; ++ix, ++j) {
            x[j] = ix;
            y[j] = iy;
        }
    }
    _impl.reset(new Impl(x, y, model->makeEllipseVector(), model->getBasisVector(), sigma));
    _data = ndarray::flatten<1>(ndarray::copy(image));
    _data.deep() /= sigma;
    _weights = ndarray::allocate(_data.getShape());
    _weights.deep() = 1.0;
}

void MultiShapeletPsfLikelihood::computeModelMatrix(
    ndarray::Array<Pixel,2,-1> const & modelMatrix,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    bool doApplyWeights
) const {
    return _impl->computeModelMatrix(modelMatrix, nonlinear, _fixed, *getModel());
}

MultiShapeletPsfLikelihood::~MultiShapeletPsfLikelihood() {}

}}} // namespace lsst::meas::modelfit
