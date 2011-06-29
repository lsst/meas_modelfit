// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include <Eigen/Cholesky>
#include <fstream>
#include "boost/serialization/binary_object.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"


namespace afwShapelets = lsst::afw::math::shapelets;

namespace lsst { namespace meas { namespace multifit {

namespace {

class ConvolvedCompoundShapeletModelBasis : public ModelBasis {
public:

    struct Element {

        ConvolvedShapeletModelBasis::Ptr component;
        ndarray::EigenView<Pixel,2,1> forward;

        Element(
            ConvolvedShapeletModelBasis::Ptr const & component_, 
            ndarray::EigenView<Pixel,2,1> const & forward_
        ) :
            component(component_), forward(forward_)
        {}

        Element & operator=(Element const & other) {
            if (&other != this) {
                component = other.component;
                forward.setArray(other.forward.getArray());
            }
            return *this;
        }
    };

    virtual Eigen::MatrixXd computeInnerProductMatrix(
        lsst::afw::geom::ellipses::BaseCore const & ellipse
    ) const {
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(getSize(), getSize());
        for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
            afw::geom::Ellipse frontEllipse1(ellipse);
            frontEllipse1.scale(i->component->getScale());
            Eigen::MatrixXd p1 = i->forward.transpose() 
                * ndarray::viewAsTransposedEigen(i->component->getConvolution()->evaluate(frontEllipse1));
            double r1 = frontEllipse1.getCore().getDeterminantRadius();
            for (ElementVector::const_iterator j = _elements.begin(); j != _elements.end(); ++j) {
                afw::geom::Ellipse frontEllipse2(ellipse);
                frontEllipse2.scale(j->component->getScale());
                Eigen::MatrixXd p2 = 
                    ndarray::viewAsEigen(j->component->getConvolution()->evaluate(frontEllipse2))
                    * j->forward;
                double r2 = frontEllipse2.getCore().getDeterminantRadius();
                Eigen::MatrixXd m = afwShapelets::HermiteEvaluator::computeInnerProductMatrix(
                    i->component->getConvolution()->getRowOrder(), 
                    j->component->getConvolution()->getRowOrder(), 
                    1.0 / r1, 1.0 / r2
                );
                m /= i->component->getScale() * j->component->getScale();
                result += p1 * m * p2;
            }
        }
        return result;
    }

    typedef std::vector<Element> ElementVector;

    ConvolvedCompoundShapeletModelBasis(int size, ElementVector & elements) : ModelBasis(size), _elements() {
        elements.swap(_elements);
    }

protected:

    virtual void _integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot integrate convolved basis."
        );      
    }

    virtual void _evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const {
        matrix.deep() = 0.0;
        for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
            ndarray::Array<Pixel,2,2> front =
                ndarray::allocate(footprint->getArea(), i->component->getSize());
            i->component->evaluate(front, footprint, ellipse);
            ndarray::viewAsEigen(matrix) += ndarray::viewAsEigen(front) * i->forward;
        }
    }

private:
    ElementVector _elements;
};


} // anonymous

namespace detail {

class CompoundShapeletImpl {
public:

    typedef boost::shared_ptr<CompoundShapeletImpl> Ptr;
    typedef std::vector<ShapeletModelBasis::Ptr> ComponentVector;

    CompoundShapeletImpl(
        ComponentVector const & components,
        ndarray::Array<Pixel,2,1> const & forward,
        ndarray::Array<Pixel,2,1> const & reverse
    );
    
    explicit CompoundShapeletImpl(ComponentVector const & components);

    CompoundShapeletImpl(CompoundShapeletImpl const & other);
    
    ComponentVector extractComponents() const;

    lsst::ndarray::Array<Pixel,2,1> getForward() const { return _forward; }
    lsst::ndarray::Array<Pixel,2,1> getReverse() const { return _reverse; }

    int getSize() const { return _forward.getSize<1>(); }

    ModelBasis::Ptr convolve(lsst::afw::math::shapelets::ShapeletFunction const & psf) const;
    ModelBasis::Ptr convolve(lsst::afw::math::shapelets::MultiShapeletFunction const & psf) const;
    
    void integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const;

    void evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const;

    void evaluateRadialProfile(
        lsst::ndarray::Array<Pixel,2,1> const & profile,
        lsst::ndarray::Array<Pixel const,1,1> const & radii
    ) const;

    Eigen::MatrixXd computeInnerProductMatrix(
        lsst::afw::geom::ellipses::BaseCore const & ellipse
    ) const;

    Eigen::MatrixXd computeInnerProductMatrix() const;

    void normalize();

    void orthogonalize();

    void slice(int start, int stop);

    void setMapping(
        lsst::ndarray::Array<Pixel const,2,1> const & forward,
        lsst::ndarray::Array<Pixel const,2,1> const & reverse
    );

    static Ptr load(std::string const & filename);
    void save(std::string const & filename);

protected:

    typedef ndarray::EigenView<Pixel,2,1> Matrix;
    typedef ndarray::TransposedEigenView<Pixel,2,1> MatrixT;

    struct Element {
        ShapeletModelBasis::Ptr component;
        Matrix forward;
        MatrixT reverse;

        Element(
            ShapeletModelBasis::Ptr const & component_, 
            ndarray::Array<Pixel,2,1> const & fullForward,
            ndarray::Array<Pixel,2,1> const & fullReverse,
            int offset
        ) :
            component(component_),
            forward(fullForward[ndarray::view(offset, offset + component->getSize())()]),
            reverse(fullReverse[ndarray::view(offset, offset + component->getSize())()])
        {}

        Element & operator=(Element const & other) {
            if (&other != this) {
                component = other.component;
                forward.setArray(other.forward.getArray());
                reverse.setArray(other.reverse.getArray());
            }
            return *this;
        }
    };

    typedef std::vector<Element> ElementVector;

    void _fillElements(ComponentVector const & components);
    void _resetElements();

    static int _computeSize(ComponentVector const & components);

    static ndarray::Array<Pixel,2,2> _makeIdentity(int size);

    ElementVector _elements;
    ndarray::Array<Pixel,2,1> _forward;
    ndarray::Array<Pixel,2,1> _reverse;
};

CompoundShapeletImpl::CompoundShapeletImpl(
    ComponentVector const & components, 
    ndarray::Array<Pixel,2,1> const & fullForward,
    ndarray::Array<Pixel,2,1> const & fullReverse
) :
    _elements(),
    _forward(fullForward),
    _reverse(fullReverse)
{
    checkSize(
        _computeSize(components), _forward.getSize<0>(),
        "Aggregate component shapelet basis size (%d) does not match forward mapping rows (%d)."
    );
    checkSize(
        _reverse.getSize<0>(), _forward.getSize<0>(),
        "Reverse mapping rows (%d) does not match forward mapping rows (%d)."
    );
    checkSize(
        _reverse.getSize<1>(), _forward.getSize<1>(),
        "Reverse mapping columns (%d) does not match forward mapping columns (%d)."
    );
    _fillElements(components);
}

CompoundShapeletImpl::CompoundShapeletImpl(ComponentVector const & components) :
    _elements(),
    _forward(_makeIdentity(_computeSize(components))),
    _reverse(ndarray::copy(_forward))
{
    _fillElements(components);
}

CompoundShapeletImpl::CompoundShapeletImpl(CompoundShapeletImpl const & other) :
    _elements(other._elements),
    _forward(ndarray::copy(other._forward)),
    _reverse(ndarray::copy(other._reverse))
{
    _resetElements();
}

CompoundShapeletImpl::ComponentVector 
CompoundShapeletImpl::extractComponents() const {
    ComponentVector result;
    result.reserve(_elements.size());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        result.push_back(i->component);
    }
    return result;
}

Eigen::MatrixXd CompoundShapeletImpl::computeInnerProductMatrix() const {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_forward.getSize<1>(), _forward.getSize<1>());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        for (ElementVector::const_iterator j = _elements.begin(); j != _elements.end(); ++j) {
            Eigen::MatrixXd m = afwShapelets::HermiteEvaluator::computeInnerProductMatrix(
                i->component->getOrder(), j->component->getOrder(),
                1.0 / i->component->getScale(), 1.0 / j->component->getScale()
            );
            m /= (i->component->getScale() * j->component->getScale());
            result += i->forward.transpose() * m * j->forward;
        }
    }
    return result;
}

Eigen::MatrixXd CompoundShapeletImpl::computeInnerProductMatrix(
    lsst::afw::geom::ellipses::BaseCore const & ellipse
) const {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_forward.getSize<1>(), _forward.getSize<1>());
    double r = ellipse.getDeterminantRadius();
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        for (ElementVector::const_iterator j = _elements.begin(); j != _elements.end(); ++j) {
            Eigen::MatrixXd m = afwShapelets::HermiteEvaluator::computeInnerProductMatrix(
                i->component->getOrder(), j->component->getOrder(),
                1.0 / (r * i->component->getScale()), 1.0 / (r * j->component->getScale())
            );
            result += i->forward.transpose() * m * j->forward;
        }
    }
    return result;
}

ModelBasis::Ptr CompoundShapeletImpl::convolve(
    afwShapelets::ShapeletFunction const & psf
) const {
    ConvolvedCompoundShapeletModelBasis::ElementVector convolvedElements;
    convolvedElements.reserve(_elements.size());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        convolvedElements.push_back(
            ConvolvedCompoundShapeletModelBasis::Element(i->component->convolve(psf), i->forward)
        );
    }
    return ModelBasis::Ptr(
        new ConvolvedCompoundShapeletModelBasis(this->getSize(), boost::ref(convolvedElements))
    );
}

ModelBasis::Ptr CompoundShapeletImpl::convolve(
    afwShapelets::MultiShapeletFunction const & psf
) const {
    ConvolvedCompoundShapeletModelBasis::ElementVector convolvedElements;
    convolvedElements.reserve(_elements.size() * psf.getElements().size());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        for (LocalPsf::MultiShapelet::ElementList::const_iterator j = psf.getElements().begin(); 
             j != psf.getElements().end();
             ++j
        ) {
            convolvedElements.push_back(
                ConvolvedCompoundShapeletModelBasis::Element(i->component->convolve(*j), i->forward)
            );
        }
    }
    return ModelBasis::Ptr(
        new ConvolvedCompoundShapeletModelBasis(this->getSize(), boost::ref(convolvedElements))
    );
}

void CompoundShapeletImpl::evaluate(
    ndarray::Array<Pixel, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
    afw::geom::Ellipse const & ellipse
) const {
    matrix.deep() = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        ndarray::Array<Pixel,2,2> front(ndarray::allocate(footprint->getArea(), i->component->getSize()));
        i->component->evaluate(front, footprint, ellipse);
        ndarray::viewAsEigen(matrix) += ndarray::viewAsEigen(front) * i->forward;
    }
}

void CompoundShapeletImpl::integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    vector.deep() = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        ndarray::Array<Pixel,1,1> front(ndarray::allocate(i->component->getSize()));
        i->component->integrate(front);
        ndarray::viewAsTransposedEigen(vector) += ndarray::viewAsTransposedEigen(front) * i->forward;
    }
}

void CompoundShapeletImpl::evaluateRadialProfile(
    lsst::ndarray::Array<Pixel,2,1> const & profile,
    lsst::ndarray::Array<Pixel const,1,1> const & radii
) const {
    profile.deep() = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        ndarray::Array<Pixel,2,2> front(ndarray::allocate(radii.getSize<0>(), i->component->getSize()));
        i->component->evaluateRadialProfile(front, radii);
        ndarray::viewAsEigen(profile) += ndarray::viewAsEigen(front) * i->forward;
    }
}

void CompoundShapeletImpl::normalize() {
    Eigen::MatrixXd v = computeInnerProductMatrix();
    for (int n = 0; n < getSize(); ++n) {
        double s = std::sqrt(v(n, n));
        ndarray::viewAsEigen(_forward).col(n) /= s;
        ndarray::viewAsEigen(_reverse).col(n) *= s;
    }
}

void CompoundShapeletImpl::orthogonalize() {
    Eigen::MatrixXd v = computeInnerProductMatrix();
    Eigen::LLT<Eigen::MatrixXd> cholesky(v);
    Eigen::MatrixXd m = Eigen::MatrixXd::Identity(v.rows(), v.cols());
    cholesky.matrixL().transpose().solveTriangularInPlace(m);
    ndarray::viewAsEigen(_forward) = ndarray::viewAsEigen(_forward) * m.part<Eigen::LowerTriangular>();
    ndarray::viewAsTransposedEigen(_reverse) 
        = cholesky.matrixL().transpose() * ndarray::viewAsEigen(_reverse);
}

void CompoundShapeletImpl::slice(int start, int stop) {
    _forward = _forward[ndarray::view()(start, stop)];
    _reverse = _reverse[ndarray::view()(start, stop)];
    _resetElements();
}

void CompoundShapeletImpl::setMapping(
    ndarray::Array<Pixel const,2,1> const & forward,
    ndarray::Array<Pixel const,2,1> const & reverse
) {
    _forward = ndarray::copy(forward);
    _reverse = ndarray::copy(reverse);
    _resetElements();
}

void CompoundShapeletImpl::_fillElements(ComponentVector const & components) {
    _elements.reserve(components.size());
    ComponentVector::const_iterator i = components.begin();
    for (int offset = 0; i != components.end(); ++i) {
        _elements.push_back(Element(*i, _forward, _reverse, offset));
        offset += (**i).getSize();
    }
}

void CompoundShapeletImpl::_resetElements() {
    ElementVector new_elements;
    new_elements.reserve(_elements.size());
    ElementVector::const_iterator i = _elements.begin();
    for (int offset = 0; i != _elements.end(); ++i) {
        new_elements.push_back(Element(i->component, _forward, _reverse, offset));
        offset += i->component->getSize();
    }
    _elements.swap(new_elements);
}

int CompoundShapeletImpl::_computeSize(ComponentVector const & components) {
    int size = 0;
    for (ComponentVector::const_iterator i = components.begin(); i != components.end(); ++i) {
        size += (**i).getSize();
    }
    return size;
}

ndarray::Array<Pixel,2,2> CompoundShapeletImpl::_makeIdentity(int size) {
    ndarray::Array<Pixel,2,2> result(ndarray::allocate(size, size));
    ndarray::viewAsEigen(result).setIdentity();
    return result;
}

CompoundShapeletImpl::Ptr CompoundShapeletImpl::load(
    std::string const & filename
) {
    std::ifstream ifs(filename.c_str());
    boost::archive::text_iarchive ar(ifs);
    
    int nComponents;
    ar >> nComponents;
    ComponentVector components(nComponents);        
    int order;
    double scale;
    for(int i =0; i < nComponents; ++i) {
        ar >> order;
        ar >> scale;
        components[i] = multifit::ShapeletModelBasis::make(order, scale);            
    }
    int width, height;
    ar >> height;
    ar >> width;
    int size = width*height;
    ndarray::Array<Pixel, 2, 2> forward = ndarray::allocate(
        ndarray::makeVector(height, width)
    );
    ar >> boost::serialization::make_array(forward.getData(), size);
    ndarray::Array<Pixel, 2, 2> reverse = ndarray::allocate(
        ndarray::makeVector(height, width)
    );
    ar >> boost::serialization::make_array(reverse.getData(), size);
    return CompoundShapeletImpl::Ptr(new CompoundShapeletImpl(components, forward, reverse));
}

void CompoundShapeletImpl::save(std::string const & filename) {
    std::ofstream ofs(filename.c_str());
    boost::archive::text_oarchive ar(ofs);

    int nElement = _elements.size();
    ar << nElement;
    for(int i=0; i < nElement; ++i) {
        int order = _elements[i].component->getOrder();
        double scale = _elements[i].component->getScale();
        ar << order;
        ar << scale;
    }
    int height = _forward.getSize<0>(), width = _forward.getSize<1>();
    int size = width*height;
    ndarray::Array<Pixel, 2, 2> forward = ndarray::copy(_forward);
    ndarray::Array<Pixel, 2, 2> reverse = ndarray::copy(_reverse);
    ar << height;
    ar << width;
    ar << boost::serialization::make_array(forward.getData(), size);
    ar << boost::serialization::make_array(reverse.getData(), size);
}

} //namespace detail

CompoundShapeletModelBasis::ComponentVector CompoundShapeletModelBasis::extractComponents() const {
    return _impl->extractComponents();
}

ndarray::Array<Pixel const,2,1> CompoundShapeletModelBasis::getForward() const {
    return _impl->getForward();
}

ndarray::Array<Pixel const,2,1> CompoundShapeletModelBasis::getReverse() const {
    return _impl->getReverse();
}

ModelBasis::Ptr CompoundShapeletModelBasis::convolve(
    LocalPsf::ConstPtr const & psf
) const {
    if (psf->hasNativeShapelet()) {
        afwShapelets::MultiShapeletFunction s = psf->getNativeShapelet(afwShapelets::HERMITE);
        s.shiftInPlace(-afw::geom::Extent2D(psf->getPoint()));
        return convolve(s);
    } else {
        afwShapelets::ShapeletFunction s = 
            psf->computeShapelet(afwShapelets::HERMITE, ShapeletModelBasis::getPsfShapeletOrder());
        s.getEllipse().getCenter() -= afw::geom::Extent2D(psf->getPoint());
        return convolve(s);
    }
}

ModelBasis::Ptr CompoundShapeletModelBasis::convolve(
    afwShapelets::ShapeletFunction const & psf
) const {
    return _impl->convolve(psf);
}

ModelBasis::Ptr CompoundShapeletModelBasis::convolve(
    afwShapelets::MultiShapeletFunction const & psf
) const {
    return _impl->convolve(psf);
}

Eigen::MatrixXd CompoundShapeletModelBasis::computeInnerProductMatrix() const {
    return _impl->computeInnerProductMatrix();
}

Eigen::MatrixXd CompoundShapeletModelBasis::computeInnerProductMatrix(
    lsst::afw::geom::ellipses::BaseCore const & ellipse
) const {
    return _impl->computeInnerProductMatrix(ellipse);
}

void CompoundShapeletModelBasis::_evaluate(
    ndarray::Array<Pixel, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
    afw::geom::Ellipse const & ellipse
) const {
    _impl->evaluate(matrix, footprint, ellipse);
}

void CompoundShapeletModelBasis::_integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    _impl->integrate(vector);
}

void CompoundShapeletModelBasis::_evaluateRadialProfile(
    lsst::ndarray::Array<Pixel,2,1> const & profile,
    lsst::ndarray::Array<Pixel const,1,1> const & radii
) const {
    _impl->evaluateRadialProfile(profile, radii);
}

CompoundShapeletModelBasis::CompoundShapeletModelBasis(
    detail::CompoundShapeletImpl::Ptr const & impl
) : ModelBasis(impl->getSize()), _impl(impl)
{}

CompoundShapeletModelBasis::~CompoundShapeletModelBasis() {}

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components
) : _impl(new detail::CompoundShapeletImpl(components))
{}

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components,
    ndarray::Array<Pixel const,2,1> const & forward,
    ndarray::Array<Pixel const,2,1> const & reverse
) : _impl(new detail::CompoundShapeletImpl(components, ndarray::copy(forward), ndarray::copy(reverse)))
{}

CompoundShapeletBuilder::ComponentVector CompoundShapeletBuilder::extractComponents() const {
    return _impl->extractComponents();
}

ndarray::Array<Pixel const,2,1> CompoundShapeletBuilder::getForward() const {
    return _impl->getForward();
}

ndarray::Array<Pixel const,2,1> CompoundShapeletBuilder::getReverse() const {
    return _impl->getReverse();
}

int CompoundShapeletBuilder::getSize() const {
    return _impl->getSize();
}

Eigen::MatrixXd CompoundShapeletBuilder::computeInnerProductMatrix() const {
    return _impl->computeInnerProductMatrix();
}

void CompoundShapeletBuilder::edit() {
    if (!_impl.unique()) {
        detail::CompoundShapeletImpl::Ptr impl(new detail::CompoundShapeletImpl(*_impl));
        _impl.swap(impl);
    }
}

void CompoundShapeletBuilder::normalize() {
    edit();
    _impl->normalize();
}

void CompoundShapeletBuilder::orthogonalize() {
    edit();
    _impl->orthogonalize();    
}

void CompoundShapeletBuilder::slice(int start, int stop) {
    edit();
    _impl->slice(start, stop);
}

void CompoundShapeletBuilder::setMapping(
    ndarray::Array<Pixel const,2,1> const & forward,
    ndarray::Array<Pixel const,2,1> const & reverse
) {
    detail::checkSize(
        getForward().getSize<0>(), forward.getSize<0>(),
        "Aggregate component shapelet basis size (%d) does not match forward mapping rows (%d)."
    );
    detail::checkSize(
        reverse.getSize<0>(), forward.getSize<0>(),
        "Reverse mapping rows (%d) does not match forward mapping rows (%d)."
    );
    detail::checkSize(
       reverse.getSize<1>(), forward.getSize<1>(),
        "Reverse mapping columns (%d) does not match forward mapping columns (%d)."
    );
    edit();
    _impl->setMapping(forward, reverse);
}

CompoundShapeletModelBasis::Ptr CompoundShapeletBuilder::build() const {
    return CompoundShapeletModelBasis::Ptr(new CompoundShapeletModelBasis(_impl));
}

CompoundShapeletModelBasis::Ptr CompoundShapeletModelBasis::load(
    std::string const & filename
) {
    return CompoundShapeletModelBasis::Ptr(
        new CompoundShapeletModelBasis(
            detail::CompoundShapeletImpl::load(filename)
        )
    );
}

void CompoundShapeletModelBasis::save(std::string const & filename) {
    _impl->save(filename);
}

CompoundShapeletBuilder::~CompoundShapeletBuilder() {}

}}} //end namespace lsst::meas::multifit
