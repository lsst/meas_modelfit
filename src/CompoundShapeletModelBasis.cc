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
        ndarray::EigenView<const Pixel,2,1> forward;

        Element(
            ConvolvedShapeletModelBasis::Ptr const & component_, 
            ndarray::EigenView<const Pixel,2,1> const & forward_
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
        return Eigen::MatrixXd::Identity(getSize(), getSize());
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

CompoundShapeletBase::ComponentVector 
CompoundShapeletBase::extractComponents() const {
    ComponentVector result;
    result.reserve(_elements.size());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        result.push_back(i->component);
    }
    return result;
}

Eigen::MatrixXd CompoundShapeletBase::computeInnerProductMatrix() const {
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

CompoundShapeletBase::Element::Element(
    ShapeletModelBasis::Ptr const & component_, 
    ndarray::Array<const Pixel,2,1> const & fullForward,
    ndarray::Array<const Pixel,2,1> const & fullReverse,
    int offset
) : 
    component(component_),
    forward(fullForward[ndarray::view(offset, offset + component->getSize())()]),
    reverse(fullReverse[ndarray::view(offset, offset + component->getSize())()])
{}

CompoundShapeletBase::Element & 
CompoundShapeletBase::Element::operator=(Element const & other) {
    if (&other != this) {
        component = other.component;
        forward.setArray(other.forward.getArray());
        reverse.setArray(other.reverse.getArray());
    }
    return *this;
}

CompoundShapeletBase::CompoundShapeletBase(ComponentVector const & components) :
    _elements(),
    _forward(_makeIdentity(_computeSize(components))),
    _reverse(ndarray::copy(_forward))
{
    _fillElements(components);
}

CompoundShapeletBase::CompoundShapeletBase(
    ComponentVector const & components, 
    ndarray::Array<const Pixel,2,1> const & fullForward,
    ndarray::Array<const Pixel,2,1> const & fullReverse
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

void CompoundShapeletBase::_fillElements(ComponentVector const & components) {
    _elements.reserve(components.size());
    ComponentVector::const_iterator i = components.begin();
    for (int offset = 0; i != components.end(); ++i) {
        _elements.push_back(Element(*i, _forward, _reverse, offset));
        offset += (**i).getSize();
    }
}

void CompoundShapeletBase::_resetElements() {
    ElementVector new_elements;
    new_elements.reserve(_elements.size());
    ElementVector::const_iterator i = _elements.begin();
    for (int offset = 0; i != _elements.end(); ++i) {
        new_elements.push_back(Element(i->component, _forward, _reverse, offset));
        offset += i->component->getSize();
    }
    _elements.swap(new_elements);
}

int CompoundShapeletBase::_computeSize(ComponentVector const & components) {
    int size = 0;
    for (ComponentVector::const_iterator i = components.begin(); i != components.end(); ++i) {
        size += (**i).getSize();
    }
    return size;
}

ndarray::Array<Pixel,2,2> CompoundShapeletBase::_makeIdentity(int size) {
    ndarray::Array<Pixel,2,2> result(ndarray::allocate(size, size));
    ndarray::viewAsEigen(result).setIdentity();
    return result;
}

} //namespace detail


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
    ConvolvedCompoundShapeletModelBasis::ElementVector convolvedElements;
    convolvedElements.reserve(_elements.size());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        convolvedElements.push_back(
            ConvolvedCompoundShapeletModelBasis::Element(i->component->convolve(psf), i->forward)
        );
    }
    return boost::make_shared<ConvolvedCompoundShapeletModelBasis>(
        this->getSize(), boost::ref(convolvedElements)
    );
}

ModelBasis::Ptr CompoundShapeletModelBasis::convolve(
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
    return boost::make_shared<ConvolvedCompoundShapeletModelBasis>(
        this->getSize(), boost::ref(convolvedElements)
    );
}

void CompoundShapeletModelBasis::_evaluate(
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

void CompoundShapeletModelBasis::_integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    vector.deep() = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        ndarray::Array<Pixel,1,1> front(ndarray::allocate(i->component->getSize()));
        i->component->integrate(front);
        ndarray::viewAsTransposedEigen(vector) += ndarray::viewAsTransposedEigen(front) * i->forward;
    }
}

void CompoundShapeletModelBasis::_evaluateRadialProfile(
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

CompoundShapeletModelBasis::CompoundShapeletModelBasis(
    CompoundShapeletBuilder const & builder
) : ModelBasis(builder.getSize()),
    detail::CompoundShapeletBase(builder)
{
    if (builder._constraintMatrix.getSize<0>() > 0) {
        attachConstraint(builder._constraintMatrix, builder._constraintVector);
    }
}

Eigen::MatrixXd CompoundShapeletModelBasis::computeInnerProductMatrix(
    lsst::afw::geom::ellipses::BaseCore const & ellipse
) const {
    return Eigen::MatrixXd::Identity(getSize(), getSize());
}

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components
) :
    detail::CompoundShapeletBase(components)
{}

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components,
    ndarray::Array<Pixel const,2,1> const & forward,
    ndarray::Array<Pixel const,2,1> const & reverse
) : 
    detail::CompoundShapeletBase(components, forward, reverse) 
{}

void CompoundShapeletBuilder::orthogonalize() {
    Eigen::MatrixXd v = computeInnerProductMatrix();
    Eigen::LLT<Eigen::MatrixXd> cholesky(v);
    Eigen::MatrixXd m = Eigen::MatrixXd::Identity(v.rows(), v.cols());
    cholesky.matrixL().transpose().solveTriangularInPlace(m);
    Matrix newForward(ndarray::allocate(_forward.getShape()));
    MatrixT newReverse(ndarray::allocate(_reverse.getShape()));
    newForward = ndarray::viewAsEigen(_forward) * m;
    newReverse = cholesky.matrixL().transpose() * ndarray::viewAsEigen(_reverse);
    if (_constraintMatrix.getSize<0>() > 0) {
        ndarray::Array<Pixel,2,2> newConstraintMatrix = ndarray::allocate(_constraintMatrix.getShape());
        ndarray::viewAsEigen(newConstraintMatrix) = ndarray::viewAsEigen(_constraintMatrix) * m;
    }
    _forward = newForward.getArray();
    _reverse = newReverse.getArray();
    _resetElements();
}

void CompoundShapeletBuilder::slice(int start, int stop) {
    _forward = _forward[ndarray::view()(start, stop)];
    _reverse = _reverse[ndarray::view()(start, stop)];
    _constraintMatrix = _constraintMatrix[ndarray::view()(start, stop)];
    _resetElements();

}

void CompoundShapeletBuilder::setMapping(
    ndarray::Array<Pixel const,2,1> const & forward,
    ndarray::Array<Pixel const,2,1> const & reverse
) {
    detail::checkSize(
        _forward.getSize<0>(), forward.getSize<0>(),
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
    _forward = ndarray::copy(forward);
    _reverse = ndarray::copy(reverse);
    _resetElements();
}

void CompoundShapeletBuilder::setConstraint(
    ndarray::Array<Pixel const,2,1> const & matrix,
    ndarray::Array<Pixel const,1,1> const & vector
) {
    detail::checkSize(
        matrix.getSize<0>(), vector.getSize<0>(),
        "Number of constraints in matrix (%d) do not match number of constraints in vector (%d)."
    );
    detail::checkSize(
        matrix.getSize<1>(), getSize(),
        "Incorrect number of columns (%d) in constraint matrix (expected %d)."
    );
    _constraintMatrix = ndarray::copy(matrix);
    _constraintVector = ndarray::copy(vector);
}

CompoundShapeletModelBasis::Ptr CompoundShapeletBuilder::build() const {
    return boost::make_shared<CompoundShapeletModelBasis>(*this);
}

CompoundShapeletModelBasis::Ptr CompoundShapeletModelBasis::load(
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
    int constraintSize;
    ar >> constraintSize;
    CompoundShapeletBuilder builder(components, forward, reverse);
    if (constraintSize > 0) {
        ndarray::Array<Pixel,2,2> cMatrix = ndarray::allocate(constraintSize, width);
        ndarray::Array<Pixel,1,1> cVector = ndarray::allocate(constraintSize);
        ar >> boost::serialization::make_array(cMatrix.getData(), cMatrix.getNumElements());
        ar >> boost::serialization::make_array(cVector.getData(), cVector.getNumElements());
        builder.setConstraint(cMatrix, cVector);
    }
    return builder.build();
}

void CompoundShapeletModelBasis::save(std::string const & filename) {
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
    int constraintSize = getConstraintSize();
    ar << constraintSize;
    if (constraintSize > 0) {
        ndarray::Array<Pixel,2,2> cMatrix = ndarray::copy(getConstraintMatrix());
        ndarray::Array<Pixel,1,1> cVector = ndarray::copy(getConstraintVector());
        ar << boost::serialization::make_array(cMatrix.getData(), cMatrix.getNumElements());
        ar << boost::serialization::make_array(cVector.getData(), cVector.getNumElements());
    }
}

}}} //end namespace lsst::meas::multifit
