#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/ndarray/eigen.h"
#include <limits>

namespace {

template <typename PixelT>
lsst::meas::multifit::definition::Frame makeFrame(
    PTR(lsst::afw::image::Exposure<PixelT>) const & exposure,
    lsst::meas::multifit::Footprint::Ptr const & fp
) {
    //TODO: determine what mask planes to consider when masking the footprint
    typedef lsst::meas::multifit::Pixel Pixel;
    typedef lsst::afw::image::MaskedImage<PixelT> MaskedImage;
    typename MaskedImage::Mask::Pixel bitmask=~0x0;
    lsst::meas::multifit::Footprint::Ptr maskedFp = 
        lsst::afw::detection::footprintAndMask(
            fp, exposure->getMaskedImage().getMask(), bitmask
        );

    //grab portion of exposure's MaskedImage that matches fp
    lsst::afw::geom::BoxI bbox = maskedFp->getBBox();
    MaskedImage mi(
        exposure->getMaskedImage(), 
        bbox,
        lsst::afw::image::PARENT
    );

    //construct flatten version
    lsst::ndarray::Array<Pixel, 1, 1> data = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(maskedFp->getArea())
    );
    lsst::ndarray::Array<Pixel, 1, 1> variance = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(maskedFp->getArea())
    );
    lsst::afw::detection::flattenArray(
        *maskedFp, mi.getImage()->getArray(), data
    );
    lsst::afw::detection::flattenArray(
        *maskedFp, mi.getVariance()->getArray(), variance
    );
    lsst::ndarray::EigenView<Pixel, 1, 1> weights = 
        lsst::ndarray::viewAsEigen(variance).cwise().sqrt().inverse();
    lsst::meas::multifit::definition::Frame frame(
        0, 
        exposure->getFilter().getId(), 
        exposure->getWcs(), 
        exposure->getPsf(), 
        data, 
        weights.getArray()
    );
    return frame;
}

} //end annonymous namespace


namespace lsst { namespace meas { namespace multifit {

Definition Evaluator::makeDefinition() const {
    return _grid->makeDefinition();
}

Definition Evaluator::makeDefinition(
    ndarray::Array<double const,1,1> const & parameters
) const {
    detail::checkSize(
        parameters.getSize<0>(), getParameterSize(),
        "Parameter vector size (%s) is incorrect; expected (%s)."
    );
    return _grid->makeDefinition(parameters.getData());
}

Evaluator::Ptr Evaluator::make(Definition const & definition) {
    boost::shared_ptr<Grid> grid = boost::make_shared<Grid>(definition);
    return boost::make_shared<Evaluator>(grid);
}


template<typename PixelT>
Evaluator::Ptr Evaluator::make(
    PTR(afw::image::Exposure<PixelT>) const & exposure,
    Footprint::Ptr const & fp,
    afw::geom::Point2D const & position,
    bool isVariable, bool fixPosition
) {    
    //make a point source evaluator   
    Definition psDefinition;
    psDefinition.frames.insert(::makeFrame<PixelT>(exposure, fp));
    definition::Object object = definition::Object::makeStar(0, position, isVariable);
    object.position->active = !fixPosition;
    psDefinition.objects.insert(object);
    return make(psDefinition); 
}

template<typename PixelT>
Evaluator::Ptr Evaluator::make(
    PTR(afw::image::Exposure<PixelT>) const & exposure,
    Footprint::Ptr const & fp,
    afw::geom::ellipses::Ellipse const & ellipse,
    bool fixEllipticity,
    bool fixRadius,
    bool fixPosition
) {
    //make a galaxy evaluator    
    Definition sgDefinition;
    sgDefinition.frames.insert(::makeFrame<PixelT>(exposure, fp));
    definition::Object object = definition::Object::makeGalaxy(0, ModelBasis::Ptr(), ellipse);
    object.ellipticity->active = !fixEllipticity;
    object.radius->active = !fixRadius;
    object.position->active = !fixPosition;
    sgDefinition.objects.insert(object);
    return make(sgDefinition); 
}

void Evaluator::_evaluateModelMatrix(
    ndarray::Array<double,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & param
) const {
    matrix.deep() = 0.0;
    for (
        Grid::ObjectArray::const_iterator object = _grid->objects.begin();
        object != _grid->objects.end(); ++object
    ) {
        if (object->basis) {
            lsst::afw::geom::ellipses::Ellipse ellipse = object->makeEllipse(param.getData());
            for (
                Grid::SourceArray::const_iterator source = object->sources.begin();
                source != object->sources.end(); ++source
            ) {
                int coefficientOffset = object->coefficientOffset;
                if (object->isVariable) {
                    coefficientOffset += object->coefficientCount * source->frame.frameIndex;
                } else {
                    coefficientOffset += object->coefficientCount * source->frame.filterIndex;
                }            

                lsst::ndarray::Array<double,2,1> block = matrix[
                    lsst::ndarray::view(
                        source->frame.pixelOffset, 
                        source->frame.pixelOffset + source->frame.pixelCount
                    )(
                        coefficientOffset, 
                        coefficientOffset+object->coefficientCount
                    )
                ];
                source->basis->evaluate(
                    block, 
                    source->frame.footprint, 
                    ellipse.transform(source->transform)
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
                int coefficientOffset = object->coefficientOffset;
                if (object->isVariable) {
                    coefficientOffset += object->coefficientCount * source->frame.frameIndex;
                } else {
                    coefficientOffset += object->coefficientCount * source->frame.filterIndex;
                }
                ndarray::Array<double,1,0> block = matrix[
                    lsst::ndarray::view(
                        source->frame.pixelOffset, 
                        source->frame.pixelOffset + source->frame.pixelCount
                    )(
                        coefficientOffset
                    )
                ];
                source->localPsf->evaluatePointSource(
                    *source->frame.footprint, 
                    source->transform(point),
                    block
                );
                source->frame.applyWeights(block);
            }            
        }
    }
}

void Evaluator::_evaluateModelDerivative(
    ndarray::Array<double,3,3> const & derivative,
    ndarray::Array<double const,1,1> const & param
) const {
    derivative.deep() = 0.0;
    for (
        Grid::ObjectArray::const_iterator object = _grid->objects.begin();
        object != _grid->objects.end();
        ++object
    ) {
        if (object->basis) {
            lsst::afw::geom::Ellipse ellipse = object->makeEllipse(param.getData());
            for (
                Grid::SourceArray::const_iterator source = object->sources.begin();
                source != object->sources.end();
                ++source
            ) {
                int coefficientOffset = object->coefficientOffset;
                if (object->isVariable) {
                    coefficientOffset += object->coefficientCount * source->frame.frameIndex;
                } else {
                    coefficientOffset += object->coefficientCount * source->frame.filterIndex;
                }
                ndarray::Array<double, 2, 2> fiducial = lsst::ndarray::allocate(
                    lsst::ndarray::makeVector(source->frame.pixelCount, object->coefficientCount)
                );
                source->basis->evaluate(fiducial, source->frame.footprint, ellipse.transform(source->transform));
                ndarray::Array<double,3,1> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source->frame.pixelOffset, 
                            source->frame.pixelOffset + source->frame.pixelCount
                        )(
                            coefficientOffset,
                            coefficientOffset + object->coefficientCount
                        )
                    ];
                //TODO remove magic number 5 
                //this is the number of ellipse parameters
                for (int n = 0; n < 5; ++n) {
                    std::pair<int,double> p = object->perturbEllipse(ellipse, n);
                    if (p.first < 0) continue;
                    source->basis->evaluate(
                        block[p.first],                        
                        source->frame.footprint, 
                        ellipse.transform(source->transform)
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
                int coefficientOffset = object->coefficientOffset;
                if (object->isVariable) {
                    coefficientOffset += object->coefficientCount * source->frame.frameIndex;
                } else {
                    coefficientOffset += object->coefficientCount * source->frame.filterIndex;
                }
                ndarray::Array<double, 1, 1> fiducial = source->localPsf->evaluatePointSource(
                    *source->frame.footprint,
                    source->transform(point)
                );   

                ndarray::Array<double,2,0> block = 
                    derivative[
                        ndarray::view(
                        )(
                            source->frame.pixelOffset, 
                            source->frame.pixelOffset + source->frame.pixelCount
                        )(
                            coefficientOffset
                        )
                    ];
                //TODO remove magice number 2
                //this is the number of position parameters
                for (int n = 0; n < 2; ++n) {
                    std::pair<int,double> p = object->perturbPoint(point, n);
                    if (p.first < 0) continue;
 
                    source->localPsf->evaluatePointSource(
                        *source->frame.footprint,
                        source->transform(point),
                        block[p.first]
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

void Evaluator::_writeInitialParameters(ndarray::Array<double,1,1> const & param) const {
    _grid->writeParameters(param.getData());
}

Evaluator::Evaluator(boost::shared_ptr<Grid> const & grid) :
    BaseEvaluator(grid->pixelCount, grid->coefficientCount, grid->parameterCount),
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
        if (i->weights.empty()) {
            _dataVector[
                ndarray::view(i->pixelOffset, i->pixelOffset + i->pixelCount)
            ] = i->data;
        } else {
            _dataVector[
                ndarray::view(i->pixelOffset, i->pixelOffset + i->pixelCount)
            ] = (i->data / i->weights);
        }
    }
}

}}} // namespace lsst::meas::multifit
