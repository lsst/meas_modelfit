#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

template <typename T>
struct DestroyGridElement {
    void operator()(T & p) const { p.~T(); }
};

} // anonymous

namespace grid {

class Initializer {
public:

    template <SharedElementType E>
    static void transferComponents(
        Grid const & input, Definition & output, double const * paramIter
    ) {
        typedef boost::shared_ptr< definition::SharedElement<E> > DPtr;
        typedef boost::shared_ptr< grid::SharedElement<E> > GPtr;
        typedef std::map<GPtr,DPtr> Map;
        Map unique;
        Grid::ObjectComponentArray::const_iterator gi = input.objects.begin();
        Definition::ObjectComponentSet::iterator di = output.objects.begin();
        for (; gi != input.objects.end(); ++gi, ++di) {
            GPtr gp = gi->template getComponent<E>();
            if (!gp) continue;
            std::pair<GPtr,DPtr> item(gp, DPtr());
            std::pair<typename Map::iterator,bool> r = unique.insert(item);
            if (r.second) {
                r.first->second = definition::SharedElement<E>::make(gp->getValue(), gp->isActive());
                if (paramIter != 0 && r.first->second->isActive()) {
                    detail::SharedElementTraits<E>::readParameters(
                        paramIter + gp->offset,
                        r.first->second->getValue()
                    );
                }
            }
            di->template getComponent<E>() = r.first->second;
        }
    }

    template <SharedElementType E, typename Container>
    static void transferComponents(
        Definition const & input, Grid & output, Container & container
    ) {
        typedef boost::shared_ptr< definition::SharedElement<E> > DPtr;
        typedef boost::shared_ptr< grid::SharedElement<E> > GPtr;
        typedef std::map<DPtr,GPtr> Map;
        Map unique;
        Definition::ObjectComponentSet::const_iterator di = input.objects.begin();
        Grid::ObjectComponentArray::const_iterator gi = output.objects.begin();
        for (; gi != output.objects.end(); ++gi, ++di) {
            DPtr dp = di->template getComponent<E>();
            if (!dp) continue;
            std::pair<DPtr,GPtr> item(dp, GPtr());
            std::pair<typename Map::iterator,bool> r = unique.insert(item);
            if (r.second) {
                if (dp->isActive()) {
                    r.first->second.reset(new grid::SharedElement<E>(*dp, output._parameterCount));
                    output._parameterCount += grid::SharedElement<E>::SIZE;
                    container._ptrVec.push_back(r.first->second);
                } else {
                    r.first->second.reset(new grid::SharedElement<E>(*dp, -1));
                }
            }
            GPtr const * gp;
            gi->getComponentImpl(gp);
            const_cast<GPtr &>(*gp) = r.first->second;
        }
    }

    static Definition makeDefinition(Grid const & grid, double const * paramIter) {
        Wcs::Ptr wcs;
        if (grid.getWcs()) wcs = grid.getWcs()->clone();
        Definition result(wcs);
        for (Grid::FrameArray::const_iterator i = grid.frames.begin(); i != grid.frames.end(); ++i) {
            result.frames.insert(definition::Frame(*i));
        }
        for (Grid::ObjectComponentArray::const_iterator i = grid.objects.begin(); i != grid.objects.end(); ++i) {
            result.objects.insert(definition::ObjectComponent(*i));
        }
        transferComponents<POSITION>(grid, result, paramIter);
        transferComponents<RADIUS>(grid, result, paramIter);
        transferComponents<ELLIPTICITY>(grid, result, paramIter);
        return result;
    }

    static void initializeGrid(Definition const & input, Grid & output) {
        output.objects._last = output.objects._first 
            = reinterpret_cast<grid::ObjectComponent*>(output._objectData.get());
        output.frames._last = output.frames._first
            = reinterpret_cast<grid::Frame*>(output._frameData.get());
        output.sources._last = output.sources._first
            = reinterpret_cast<grid::SourceComponent*>(output._sourceData.get());

        int constraintCount = 0;
        try {
            int frameCount = 0;
            for (
                Definition::FrameSet::const_iterator i = input.frames.begin();
                i != input.frames.end();
                ++i
            ) {
                std::pair<Grid::FilterMap::iterator,bool> r = output._filters.insert(
                    std::make_pair(i->getFilterId(), output._filterCount)
                );
                if (r.second) ++output._filterCount;
                grid::Frame * newFrame = new (output.frames._last) grid::Frame(
                    *i, output._pixelCount, r.first->second, frameCount
                );
                ++output.frames._last;
                output._pixelCount += newFrame->getPixelCount();
                ++frameCount;
            }
            for (
                Definition::ObjectComponentSet::const_iterator i = input.objects.begin(); 
                i != input.objects.end();
                ++i
            ) {
                grid::ObjectComponent * newObjectComponent = new (output.objects._last) grid::ObjectComponent(
                    *i, output._coefficientCount, frameCount, output._filterCount
                );
                ++output.objects._last;
                output._coefficientCount += newObjectComponent->getCoefficientCount();
                if (newObjectComponent->getBasis()) {
                    constraintCount += newObjectComponent->getBasis()->getConstraintSize()
                        * (newObjectComponent->getCoefficientCount() / newObjectComponent->getSourceCoefficientCount());
                } else {
                    constraintCount += 
                        (newObjectComponent->getCoefficientCount() / newObjectComponent->getSourceCoefficientCount());
                }
            }
            transferComponents<POSITION>(input, output, output.positions);
            transferComponents<RADIUS>(input, output, output.radii);
            transferComponents<ELLIPTICITY>(input, output, output.ellipticities);
            if (constraintCount) {
                output._constraintVector = ndarray::allocate(constraintCount);
                output._constraintMatrix = ndarray::allocate(constraintCount, output._coefficientCount);
                output._constraintVector.deep() = 0.0;
                output._constraintMatrix.deep() = 0.0;
            }
            int constraintOffset = 0;
            for (
                grid::ObjectComponent * i = output.objects._first; 
                i != output.objects._last;
                ++i
            ) {
                i->validate();
                i->sources._first = output.sources._last;
                for (
                    Grid::FrameArray::const_iterator j = output.frames.begin();
                    j != output.frames.end();
                    ++j
                ) {
                    new (output.sources._last) grid::SourceComponent(*j, *i, output.getWcs());
                    ++output.sources._last;
                }
                i->sources._last = output.sources._last;
                if (constraintCount) {
                    ndarray::Array<Pixel,1,1> subConstraintVector;
                    ndarray::Array<Pixel,2,1> subConstraintMatrix;
                    int nObjConstraints = 0;
                    if (i->getBasis() && i->getBasis()->getConstraintSize()) {
                        nObjConstraints = i->getBasis()->getConstraintSize();
                        subConstraintVector 
                            = ndarray::const_array_cast<Pixel>(i->getBasis()->getConstraintVector());
                        subConstraintMatrix
                            = ndarray::const_array_cast<Pixel>(i->getBasis()->getConstraintMatrix());
                    } else if (!i->getBasis()) {
                        nObjConstraints = 1;
                        subConstraintMatrix = ndarray::allocate(1,1);
                        subConstraintMatrix.deep() = 1.0;
                        subConstraintVector = ndarray::allocate(1);
                        subConstraintVector.deep() = 0.0;
                    }
                    if (nObjConstraints) {
                        int nSteps = output._filterCount;
                        if (i->isVariable()) {
                            nSteps = frameCount;
                        }
                        for (int step = 0; step < nSteps; ++step) {
                            output._constraintVector[
                                ndarray::view(constraintOffset, constraintOffset + nObjConstraints)
                            ] = subConstraintVector;
                            output._constraintMatrix[
                                ndarray::view(
                                    constraintOffset, constraintOffset + nObjConstraints
                                ) (
                                    i->getCoefficientOffset() + step * i->getSourceCoefficientCount(),
                                    i->getCoefficientOffset() + (step + 1) * i->getSourceCoefficientCount()
                                )
                            ] = subConstraintMatrix;
                            constraintOffset += nObjConstraints;
                        }
                    }
                }
            }
        } catch (...) {
            destroyGrid(output);
            throw;
        }
    }

    static void destroyGrid(Grid & g) {
        std::for_each(g.sources._first, g.sources._last, DestroyGridElement<grid::SourceComponent>());
        std::for_each(g.objects._first, g.objects._last, DestroyGridElement<grid::ObjectComponent>());
        std::for_each(g.frames._first, g.frames._last, DestroyGridElement<grid::Frame>());
    }

};

Frame::Frame(definition::Frame const & def, int pixelOffset, int filterIndex, int frameIndex) :
    detail::FrameBase(def, true), _pixelOffset(pixelOffset), _pixelCount(_footprint->getNpix()), 
    _filterIndex(filterIndex), _frameIndex(frameIndex)
{}

void Frame::applyWeights(ndarray::Array<Pixel,2,1> const & matrix) const {
    if (!_weights.empty()) {
        matrix.deep() *= _weights;
    }
}

void Frame::applyWeights(ndarray::Array<Pixel,1,0> const & vector) const {
    if (!_weights.empty()) {
        vector.deep() *= _weights;
    }
}

std::ostream & operator<<(std::ostream & os, Frame const & frame) {
    std::string filterName("undefined");
    try {
        filterName = lsst::afw::image::Filter(frame.getFilterId()).getName();
    } catch (lsst::pex::exceptions::NotFoundException &) {}
    return os << "Frame " << frame.id << " (@" << (&frame) << ") = {" << filterName 
              << ", " << frame.getFootprint()->getArea() << "pix}\n";
}

SourceComponent::SourceComponent(
    Frame const & frame_, ObjectComponent const & object_, 
    CONST_PTR(afw::image::Wcs) const & wcs
) :
    frame(frame_), object(object_), 
    _transform()
{
    if (wcs) {
        if (!frame.getWcs()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "If the definition WCS is set, all frames must have a WCS."
            );
        }
	afw::geom::Point2D point = object.getPosition()->getValue();
        _transform = frame.getWcs()->linearizeSkyToPixel(point) * wcs->linearizePixelToSky(point);
    } else {
        if (frame.getWcs()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "If the definition WCS is not set, inividual frames may not have a WCS."
            );
        }
    }
    if (frame.getPsf()) {
        _localPsf = frame.getPsf()->getLocalPsf(_transform(object.getPosition()->getValue()));
    }
    if (object.getBasis()) {
        if (frame.getPsf()) {            
            _basis = object.getBasis()->convolve(_localPsf);
        } else {
            _basis = object.getBasis();
        }
    } else {
        if (!frame.getPsf()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "All objects must have a basis if any frames do not have a PSF."
            );
        }
    }
}

Definition Grid::makeDefinition() const {
    return Initializer::makeDefinition(*this, 0);
}

Definition Grid::makeDefinition(double const * paramIter) const {
    return Initializer::makeDefinition(*this, paramIter);
}

template <typename T>
T const &
find(Array<T> const & array, ID const id) {
    typename Array<T>::const_iterator iter1 = array.begin();
    typename Array<T>::const_iterator iter2;
    int count, step;
    count = array.size();
    while (count > 0) {
        iter2 = iter1;
        step = count / 2;
        iter2 += step;
        if (iter2->id < id) {
            iter1 = ++iter2;
            count -= step + 1;
        } else {
            count = step;
        }
    }
    if (iter1->id != id) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("ObjectComponent or Frame with ID %d not found.") % id).str()
        );
    }
    return *iter1;
}

template ObjectComponent const & find(Array<ObjectComponent> const &, ID const);
template Frame const & find(Array<Frame> const &, ID const);

Grid::Grid(Definition const & definition) :
    _filterCount(0),
    _coefficientCount(0),
    _pixelCount(0),
    _parameterCount(0),
    _objectData(new char[sizeof(grid::ObjectComponent) * definition.objects.size()]),
    _frameData(new char[sizeof(grid::Frame) * definition.frames.size()]),
    _sourceData(new char[sizeof(grid::SourceComponent) * definition.frames.size() * definition.objects.size()]),
    _wcs()
{
    if (definition.getWcs()) _wcs = definition.getWcs()->clone();
    Initializer::initializeGrid(definition, *this);
}

Grid::~Grid() { Initializer::destroyGrid(*this); }

int const Grid::getFilterIndex(FilterId filterId) const {
    FilterMap::const_iterator i = _filters.find(filterId);
    if(i == _filters.end()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Filter with ID %d not found.") % filterId).str()
        );
    }
    return i->second;
}

void Grid::writeParameters(double * paramIter) const {
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        detail::SharedElementTraits<POSITION>::writeParameters(paramIter + i->offset, i->getValue());
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        detail::SharedElementTraits<RADIUS>::writeParameters(paramIter + i->offset, i->getValue());
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        detail::SharedElementTraits<ELLIPTICITY>::writeParameters(paramIter + i->offset, i->getValue());
    }
}

bool Grid::checkBounds(double const * paramIter) const {
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        if (!i->checkBounds(paramIter)) return false;
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        if (!i->checkBounds(paramIter)) return false;
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        if (!i->checkBounds(paramIter)) return false;
    }
    return true;
}

double Grid::clipToBounds(double * paramIter) const {
    double value = 0.0;
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        value += i->clipToBounds(paramIter);
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        value += i->clipToBounds(paramIter);
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        value += i->clipToBounds(paramIter);
    }
    return value;
}

}}}} // namespace lsst::meas::multifit::grid
