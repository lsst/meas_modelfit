#include "lsst/meas/multifit/Evaluator.h"

namespace lsst { namespace meas { namespace multifit {

void Evaluator::_evaluateModelMatrix(
    ndarray::Array<double,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & param
) const {
    matrix.deep() = 0.0;
    for (
        Grid::ObjectArray::const_iterator object = _grid->objects.begin();
        object != _grid->objects.end(); ++object
    ) {
        if (object->getBasis()) {
            lsst::afw::geom::ellipses::Ellipse ellipse = object->makeEllipse(param.getData());
            for (
                Grid::SourceArray::const_iterator source = object->sources.begin();
                source != object->sources.end(); ++source
            ) {
                int coefficientOffset = source->getCoefficientOffset();
                lsst::ndarray::Array<double,2,1> block = matrix[
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
                source->frame.applyWeights(block);
            }
        } else {
            afw::geom::Point2D point = object->makePoint(param.getData());
            for (
                Grid::SourceArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                ndarray::Array<double,1,0> block = matrix[
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
                source->frame.applyWeights(block);
            }            
        }
    }
}


void Evaluator::_evaluateModelMatrixDerivative(
    ndarray::Array<double,3,3> const & derivative,
    ndarray::Array<double const,2,2> const & modelMatrix,
    ndarray::Array<double const,1,1> const & param
) const {
    derivative.deep() = 0.0;
    for (
        Grid::ObjectArray::const_iterator object = _grid->objects.begin();
        object != _grid->objects.end();
        ++object
    ) {
        if (object->getBasis()) {
            lsst::afw::geom::Ellipse ellipse = object->makeEllipse(param.getData());
            for (
                Grid::SourceArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                int coefficientOffset = source->getCoefficientOffset();
                ndarray::Array<double const,2,1> fiducial = modelMatrix[
                    ndarray::view(
                        source->frame.getPixelOffset(),
                        source->frame.getPixelOffset() + source->frame.getPixelCount()
                    )(
                        coefficientOffset,
                        coefficientOffset + source->getCoefficientCount()
                    )
                ];
                ndarray::Array<double,3,1> block = 
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
            lsst::afw::geom::Point2D point = object->makePoint(param.getData());
            for (
                Grid::SourceArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                ndarray::Array<double const,1,0> fiducial = modelMatrix[
                    ndarray::view(
                        source->frame.getPixelOffset(), 
                        source->frame.getPixelOffset() + source->frame.getPixelCount()
                    )(
                        source->getCoefficientOffset()
                    )                        
                ];
                ndarray::Array<double,2,0> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source->frame.getPixelOffset(), 
                            source->frame.getPixelOffset() + source->frame.getPixelCount()
                        )(
                            source->getCoefficientOffset()
                        )
                    ];
                for (int n = 0; n < grid::PositionComponent::SIZE; ++n) {
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

Evaluator::Evaluator(Grid::Ptr const & grid) :
    BaseEvaluator(
        grid->getPixelCount(), grid->getCoefficientCount(), grid->getParameterCount(),
        -2*grid->sumLogWeights()
    ),
    _grid(grid)
{
    _initialize();
}

Evaluator::Evaluator(Evaluator const & other) : BaseEvaluator(other), _grid(other._grid) {
    _initialize();
}

void Evaluator::_initialize() {
    for (
        Grid::FrameArray::const_iterator i = _grid->frames.begin();
        i != _grid->frames.end(); ++i
    ) {
        if (i->getWeights().empty()) {
            _dataVector[
                ndarray::view(i->getPixelOffset(), i->getPixelOffset() + i->getPixelCount())
            ] = i->getData();
        } else {
            _dataVector[
                ndarray::view(i->getPixelOffset(), i->getPixelOffset() + i->getPixelCount())
            ] = (i->getData() / i->getWeights());
        }
    }
}

void Evaluator::_writeInitialParameters(ndarray::Array<double,1,1> const & param) const {
    _grid->writeParameters(param.getData());
}

}}} // namespace lsst::meas::multifit
