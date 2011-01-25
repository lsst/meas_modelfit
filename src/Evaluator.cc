#include "lsst/meas/multifit/Evaluator.h"
#include <limits>

namespace lsst { namespace meas { namespace multifit {

Definition Evaluator::makeDefinition() const {
    return _grid->makeDefinition();
}

Definition Evaluator::makeDefinition(ndarray::Array<double const,1,1> const & parameters) const {
    checkSize(
        parameters.getSize<0>(), getParameterSize(),
        "Parameter vector size (%s) is incorrect; expected (%s)."
    );
    return _grid->makeDefinition(parameters.getData());
}

Evaluator::Ptr Evaluator::make(Definition const & definition) {
    boost::shared_ptr<Grid> grid = boost::make_shared<Grid>(definition);
    return boost::make_shared<Evaluator>(grid);
}

void Evaluator::_evaluateModelMatrix(
    ndarray::Array<double,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & param
) const {
    matrix.deep() = 0.0;
    for (
        Grid::ObjectArray::const_iterator object_iter = _grid->objects.begin();
        object_iter != _grid->objects.end();
        ++object_iter
    ) {
        if (object_iter->basis) {
            lsst::afw::geom::Ellipse ellipse = object_iter->makeEllipse(param.getData());
            for (
                Grid::SourceArray::const_iterator source_iter = object_iter->sources.begin();
                source_iter != object_iter->sources.end();
                ++source_iter
            ) {
                int coefficient_offset = object_iter->coefficient_offset;
                if (object_iter->is_variable) {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.frame_index;
                } else {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.filter_index;
                }
                lsst::afw::geom::FootprintMatrix block(
                    source_iter->frame.footprint,
                    matrix[
                        ndarray::view(
                            source_iter->frame.pixel_offset, 
                            source_iter->frame.pixel_offset + source_iter->frame.pixel_count
                        )(
                            coefficient_offset,
                            coefficient_offset + object_iter->coefficient_count
                        )
                    ]
                );
                source_iter->basis->evaluate(block, ellipse.transform(source_iter->transform));
                source_iter->frame.applyWeights(block.getArray());
            }
        } else {
            lsst::afw::geom::Point2D point = object_iter->makePoint(param.getData());
            for (
                Grid::SourceArray::const_iterator source_iter = object_iter->sources.begin();
                source_iter != object_iter->sources.end();
                ++source_iter
            ) {
                int coefficient_offset = object_iter->coefficient_offset;
                if (object_iter->is_variable) {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.frame_index;
                } else {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.filter_index;
                }
                lsst::afw::detection::FootprintMatrix block(
                    source_iter->frame.footprint,
                    matrix[
                        ndarray::view(
                            source_iter->frame.pixel_offset, 
                            source_iter->frame.pixel_offset + source_iter->frame.pixel_count
                        )(
                            coefficient_offset
                        )
                    ]
                );
                source_iter->frame.psf->evaluatePointSource(block, source_iter->transform(point));
                source_iter->frame.applyWeights(block.getArray());
            }            
        }
    }
}

void Evaluator::_evaluateModelDerivative(
    ndarray::Array<double,3,3> const & derivative,
    parameters::ConstArray const & param
) const {
    derivative.deep() = 0.0;
    for (
        Grid::ObjectArray::const_iterator object_iter = _grid->objects.begin();
        object_iter != _grid->objects.end();
        ++object_iter
    ) {
        if (object_iter->basis) {
            lsst::afw::geom::Ellipse ellipse = object_iter->makeEllipse(param.getData());
            for (
                Grid::SourceArray::const_iterator source_iter = object_iter->sources.begin();
                source_iter != object_iter->sources.end();
                ++source_iter
            ) {
                int coefficient_offset = object_iter->coefficient_offset;
                if (object_iter->is_variable) {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.frame_index;
                } else {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.filter_index;
                }
                lsst::afw::detection::FootprintMatrix fiducial(
                    source_iter->frame.footprint,
                    ndarray::makeVector(object_iter->coefficient_count)
                );
                source_iter->basis->evaluate(fiducial, ellipse.transform(source_iter->transform));
                ndarray::Array<double,3,1> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source_iter->frame.pixel_offset, 
                            source_iter->frame.pixel_offset + source_iter->frame.pixel_count
                        )(
                            coefficient_offset,
                            coefficient_offset + object_iter->coefficient_count
                        )
                    ];
                for (int n = 0; n < 5; ++n) {
                    std::pair<int,double> p = object_iter->perturbEllipse(ellipse, n);
                    if (p.first < 0) continue;
                    source_iter->basis->evaluate(
                        source_iter->frame.footprint.attachArray(block[p.first].shallow()), 
                        ellipse.transform(source_iter->transform)
                    );
                    block[p.first] -= fiducial.getArray();
                    block[p.first] /= p.second;
                    object_iter->unperturbEllipse(ellipse, n, p.second);
                    source_iter->frame.applyWeights(block[p.first]);
                }
            }
        } else {
            lsst::afw::geom::Point2D point = object_iter->makePoint(param.getData());
            for (
                Grid::SourceArray::const_iterator source_iter = object_iter->sources.begin();
                source_iter != object_iter->sources.end();
                ++source_iter
            ) {
                int coefficient_offset = object_iter->coefficient_offset;
                if (object_iter->is_variable) {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.frame_index;
                } else {
                    coefficient_offset += object_iter->coefficient_count * source_iter->frame.filter_index;
                }
                lsst::afw::detection::FootprintMatrix fiducial(
                    source_iter->frame.footprint,
                    ndarray::Vector<int,0>()
                );
                source_iter->frame.psf->evaluatePointSource(fiducial, source_iter->transform(point));
                ndarray::Array<double,2,0> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source_iter->frame.pixel_offset, 
                            source_iter->frame.pixel_offset + source_iter->frame.pixel_count
                        )(
                            coefficient_offset
                        )
                    ];
                for (int n = 0; n < 2; ++n) {
                    std::pair<int,double> p = object_iter->perturbPoint(point, n);
                    if (p.first < 0) continue;
                    source_iter->frame.psf->evaluatePointSource(
                        source_iter->frame.footprint.attachArray(block[p.first].shallow()),
                        source_iter->transform(point)
                    );
                    block[p.first] -= fiducial.getArray();
                    block[p.first] /= p.second;
                    object_iter->unperturbPoint(point, n, p.second);
                    source_iter->frame.applyWeights(block[p.first]);
                }
            }
        }
    }
}

void Evaluator::_writeInitialParameters(ndarray::Array<double,1,1> const & param) const {
    _grid->writeParameters(param.getData());
}

Evaluator::Evaluator(boost::shared_ptr<Grid> const & grid) :
    BaseEvaluator(grid->pixel_count, grid->coefficient_count, grid->parameter_count),
    _grid(grid)
{
    _initialize();
}

Evaluator::Evaluator(Evaluator const & other) : BaseEvaluator(other), _grid(other._grid) {
    _initialize();
}

void Evaluator::_initialize() {
    for (
        Grid::FrameArray::const_iterator frame_iter = _grid->frames.begin();
        frame_iter != _grid->frames.end();
        ++frame_iter
    ) {
        if (frame_iter->weights.empty()) {
            _data_vector[
                ndarray::view(frame_iter->pixel_offset, frame_iter->pixel_offset + frame_iter->pixel_count)
            ] = frame_iter->data;
        } else {
            _data_vector[
                ndarray::view(frame_iter->pixel_offset, frame_iter->pixel_offset + frame_iter->pixel_count)
            ] = frame_iter->data / frame_iter->weights;
        }
    }
}

}}} // namespace lsst::meas::multifit
