#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

BaseCoefficientPrior::ConstPtr Evaluator::_evaluate(
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
                source->frame.applyWeights(block);
            }            
        }
    }
    return BaseCoefficientPrior::Ptr();
}

Evaluator::Evaluator(Grid::Ptr const & grid) :
    _grid(grid), _logPixelErrorSum(0.0),
    _dataVector(ndarray::allocate(grid->getPixelCount()))
{
    for (
        Grid::FrameArray::const_iterator i = _grid->frames.begin();
        i != _grid->frames.end(); ++i
    ) {
        _dataVector[
            ndarray::view(i->getPixelOffset(), i->getPixelOffset() + i->getPixelCount())
            ] = i->getData();
        if (!i->getWeights().empty()) {
            _dataVector[
                ndarray::view(i->getPixelOffset(), i->getPixelOffset() + i->getPixelCount())
            ] *= i->getWeights();
            _logPixelErrorSum 
                += (ndarray::viewAsEigen(i->getWeights()) * (2.0 * M_PI)).cwise().log().sum();
        }
    }
}

Evaluator::Evaluator(Evaluator const & other) :
    _grid(other._grid), _logPixelErrorSum(other._logPixelErrorSum), _dataVector(other._dataVector)
{}

void Evaluator::_writeInitialParameters(ndarray::Array<double,1,1> const & parameters) const {
    _grid->writeParameters(parameters);
}

}}} // namespace lsst::meas::multifit
