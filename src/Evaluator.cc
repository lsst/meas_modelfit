#include "lsst/meas/multifit/Evaluator.h"

namespace lsst { namespace meas { namespace multifit {

void Evaluator::_evaluateModelMatrix(
    ndarray::Array<Pixel,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & parameters
) const {
    matrix.deep() = 0.0;
    for (
        Grid::ObjectComponentArray::const_iterator object = _grid->objects.begin();
        object != _grid->objects.end(); ++object
    ) {
        if (object->getBasis()) {
            lsst::afw::geom::ellipses::Ellipse ellipse = object->makeEllipse(parameters);
            for (
                Grid::SourceComponentArray::const_iterator source = object->sources.begin();
                source != object->sources.end(); ++source
            ) {
                int coefficientOffset = source->getCoefficientOffset();
                lsst::ndarray::Array<Pixel,2,1> block = matrix[
                    lsst::ndarray::view(
                        source->frame.getPixelOffset(), 
                        source->frame.getPixelOffset() + source->frame.getPixelCount()
                    )(
                        coefficientOffset, 
                        coefficientOffset + source->getCoefficientCount()
                    )
                ];
                source->getBasis()->evaluate(
                    block, 
                    source->frame.getFootprint(), 
                    ellipse.transform(source->getTransform())
                );
                if(_usePixelWeights)
                    source->frame.applyWeights(block);
            }
        } else {
            afw::geom::Point2D point = object->makePoint(parameters);
            for (
                Grid::SourceComponentArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                ndarray::Array<Pixel,1,0> block = matrix[
                    lsst::ndarray::view(
                        source->frame.getPixelOffset(), 
                        source->frame.getPixelOffset() + source->frame.getPixelCount()
                    )(
                        source->getCoefficientOffset()
                    )
                ];
                source->getLocalPsf()->evaluatePointSource(
                    *source->frame.getFootprint(), 
                    block, 
                    source->getTransform()(point) - source->getReferencePoint()
                );
                if(_usePixelWeights)
                    source->frame.applyWeights(block);
            }            
        }
    }
}

#if 0
void Evaluator::_evaluateModelMatrixDerivative(
    ndarray::Array<Pixel,3,3> const & derivative,
    ndarray::Array<Pixel const,2,2> const & modelMatrix,
    ndarray::Array<double const,1,1> const & parameters
) const {
    derivative.deep() = 0.0;
    for (
        Grid::ObjectComponentArray::const_iterator object = _grid->objects.begin();
        object != _grid->objects.end();
        ++object
    ) {
        if (object->getBasis()) {
            lsst::afw::geom::Ellipse ellipse = object->makeEllipse(parameters);
            for (
                Grid::SourceComponentArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                int coefficientOffset = source->getCoefficientOffset();
                ndarray::Array<Pixel const,2,1> fiducial = modelMatrix[
                    ndarray::view(
                        source->frame.getPixelOffset(),
                        source->frame.getPixelOffset() + source->frame.getPixelCount()
                    )(
                        coefficientOffset,
                        coefficientOffset + source->getCoefficientCount()
                    )
                ];
                ndarray::Array<Pixel,3,1> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source->frame.getPixelOffset(), 
                            source->frame.getPixelOffset() + source->frame.getPixelCount()
                        )(
                            coefficientOffset,
                            coefficientOffset + source->getCoefficientCount()
                        )
                    ];
                //TODO remove magic number 5 
                //this is the number of ellipse parameters
                for (int n = 0; n < 5; ++n) {
                    std::pair<int,double> p = object->perturbEllipse(ellipse, n);
                    if (p.first < 0) continue;
                    source->getBasis()->evaluate(
                        block[p.first],                        
                        source->frame.getFootprint(), 
                        ellipse.transform(source->getTransform())
                    );
                    block[p.first] -= fiducial;
                    block[p.first] /= p.second;
                    object->unperturbEllipse(ellipse, n, p.second);
                    source->frame.applyWeights(block[p.first]);
                }
            }
        } else {
            lsst::afw::geom::Point2D point = object->makePoint(parameters);
            for (
                Grid::SourceComponentArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                ndarray::Array<Pixel const,1,0> fiducial = modelMatrix[
                    ndarray::view(
                        source->frame.getPixelOffset(), 
                        source->frame.getPixelOffset() + source->frame.getPixelCount()
                    )(
                        source->getCoefficientOffset()
                    )                        
                ];
                ndarray::Array<Pixel,2,0> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source->frame.getPixelOffset(), 
                            source->frame.getPixelOffset() + source->frame.getPixelCount()
                        )(
                            source->getCoefficientOffset()
                        )
                    ];
                for (int n = 0; n < grid::PositionElement::SIZE; ++n) {
                    std::pair<int,double> p = object->perturbPoint(point, n);
                    if (p.first < 0) continue;
 
                    source->getLocalPsf()->evaluatePointSource(
                        *source->frame.getFootprint(),
                        block[p.first],
                        source->getTransform()(point) - source->getReferencePoint()
                    );
                    block[p.first] -= fiducial;
                    block[p.first] /= p.second;
                    object->unperturbPoint(point, n, p.second);
                    source->frame.applyWeights(block[p.first]);
                }
            }
        }
    }
}

#endif

Evaluator::Evaluator(Grid::Ptr const & grid, bool usePixelWeights) :
    BaseEvaluator(
        grid->getPixelCount(), grid->getCoefficientCount(), grid->getParameterCount(), 
        grid->getConstraintCount()
    ),
    _usePixelWeights(usePixelWeights), _grid(grid)
{
    _constraintMatrix.deep() = 0.0;
    _constraintVector.deep() = 0.0;
    _initialize();
}

void Evaluator::_initialize() {
    for (Grid::FrameArray::iterator i = _grid->frames.begin(); i != _grid->frames.end(); ++i) {
        _dataVector[
            ndarray::view(i->getPixelOffset(), i->getPixelOffset() + i->getPixelCount())
            ] = i->getData();
        
        if (!i->getWeights().empty() && _usePixelWeights) {
            _dataVector[
                ndarray::view(i->getPixelOffset(), i->getPixelOffset() + i->getPixelCount())
            ] *= i->getWeights();
        }
    }
    for (Grid::FluxGroupArray::iterator i = _grid->groups.begin(); i != _grid->groups.end(); ++i) {
        int const nFluxInstances = i->isVariable() ? _grid->frames.size() : _grid->getFilterCount();
        for (int n = 0; n < nFluxInstances; ++n) {
            int constraintOffset = i->getConstraintOffset(n);
            int coefficientOffset = i->getCoefficientOffset(n);
            _constraintMatrix[
                ndarray::view(
                    constraintOffset, constraintOffset + i->getConstraintCount()
                )(
                    coefficientOffset, coefficientOffset + i->getSourceCoefficientCount()
                )
            ] = i->getConstraintMatrix();
            _constraintVector[
                ndarray::view(
                    constraintOffset, constraintOffset + i->getConstraintCount()
                )
            ] = i->getConstraintVector();
        }
    }
}

void Evaluator::_writeInitialParameters(ndarray::Array<double,1,1> const & parameters) const {
    _grid->writeParameters(parameters);
}

}}} // namespace lsst::meas::multifit
