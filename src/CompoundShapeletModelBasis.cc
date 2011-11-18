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
#include "lsst/meas/multifit/qp.h"
#include "Eigen/Cholesky"
#include <fstream>
#include "boost/serialization/binary_object.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

namespace afwShapelets = lsst::afw::math::shapelets;

namespace lsst { namespace meas { namespace multifit {

//================================ ConvolvedCompoundShapeletModelBasis =====================================

namespace {

class ConvolvedCompoundShapeletModelBasis : public ModelBasis {
public:

    struct Element {

        ModelBasis::Ptr component;
        ndarray::EigenView<Pixel,2,1> mapping;

        Element(ModelBasis::Ptr const & component_, ndarray::EigenView<Pixel,2,1> const & mapping_) :
            component(component_), mapping(mapping_)
        {}

        Element & operator=(Element const & other) {
            if (&other != this) {
                component = other.component;
                mapping.reset(other.mapping.shallow());
            }
            return *this;
        }
    };

    typedef std::vector<Element> ElementVector;

    ConvolvedCompoundShapeletModelBasis(int size, ElementVector & elements) : ModelBasis(size), _elements() {
        elements.swap(_elements);
    }

protected:

    virtual void _integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
        vector.deep() = 0.0;
        for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
            ndarray::Array<Pixel,1,1> front(ndarray::allocate(i->component->getSize()));
            i->component->integrate(front);
            vector.asEigen().transpose() += front.asEigen().transpose() * i->mapping;
        }
    }

    virtual void _evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        matrix.deep() = 0.0;
        for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
            ndarray::Array<Pixel,2,2> front =
                ndarray::allocate(footprint->getArea(), i->component->getSize());
            i->component->evaluate(front, footprint, ellipse);
            matrix.asEigen() += front.asEigen() * i->mapping;
        }
    }

private:
    ElementVector _elements;
};


} // anonymous

//========================================== CompoundShapeletBase =========================================

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
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_mapping.getSize<1>(), _mapping.getSize<1>());
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        for (ElementVector::const_iterator j = _elements.begin(); j != _elements.end(); ++j) {
            Eigen::MatrixXd m = afwShapelets::detail::HermiteEvaluator::computeInnerProductMatrix(
                i->component->getOrder(), j->component->getOrder(),
                1.0 / i->component->getScale(), 1.0 / j->component->getScale()
            );
            m /= (i->component->getScale() * j->component->getScale());
            result += i->mapping.transpose() * m * j->mapping;
        }
    }
    return result;
}

void CompoundShapeletBase::integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    vector.deep() = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        ndarray::Array<Pixel,1,1> front(ndarray::allocate(i->component->getSize()));
        i->component->integrate(front);
        vector.asEigen().transpose() += front.asEigen().transpose() * i->mapping;
    }
}

CompoundShapeletBase::Element::Element(
    ShapeletModelBasis::Ptr const & component_, 
    ndarray::Array<Pixel,2,1> const & fullMapping,
    int offset
) : 
    component(component_),
    mapping(fullMapping[ndarray::view(offset, offset + component->getSize())()])
{}

CompoundShapeletBase::Element & 
CompoundShapeletBase::Element::operator=(Element const & other) {
    if (&other != this) {
        component = other.component;
        mapping.reset(other.mapping.shallow());
    }
    return *this;
}

CompoundShapeletBase::CompoundShapeletBase(ComponentVector const & components) :
    _elements(),
    _mapping(_makeIdentity(_computeSize(components)))
{
    _fillElements(components);
}

CompoundShapeletBase::CompoundShapeletBase(
    ComponentVector const & components, 
    ndarray::Array<Pixel,2,1> const & fullMapping
) :
    _elements(),
    _mapping(fullMapping)
{
    checkSize(
        _computeSize(components), _mapping.getSize<0>(),
        "Aggregate component shapelet basis size (%d) does not match mapping mapping rows (%d)."
    );
    _fillElements(components);
}

void CompoundShapeletBase::_fillElements(ComponentVector const & components) {
    _elements.reserve(components.size());
    ComponentVector::const_iterator i = components.begin();
    for (int offset = 0; i != components.end(); ++i) {
        _elements.push_back(Element(*i, _mapping, offset));
        offset += (**i).getSize();
    }
}

void CompoundShapeletBase::_resetElements() {
    ElementVector new_elements;
    new_elements.reserve(_elements.size());
    ElementVector::const_iterator i = _elements.begin();
    for (int offset = 0; i != _elements.end(); ++i) {
        new_elements.push_back(Element(i->component, _mapping, offset));
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
    result.asEigen().setIdentity();
    return result;
}

} //namespace detail

//===================================== CompoundShapeletModelBasis =========================================

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
            ConvolvedCompoundShapeletModelBasis::Element(i->component->convolve(psf), i->mapping)
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
                ConvolvedCompoundShapeletModelBasis::Element(i->component->convolve(*j), i->mapping)
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
    afw::geom::ellipses::Ellipse const & ellipse
) const {
    matrix.deep() = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        ndarray::Array<Pixel,2,2> front(ndarray::allocate(footprint->getArea(), i->component->getSize()));
        i->component->evaluate(front, footprint, ellipse);
        matrix.asEigen() += front.asEigen() * i->mapping;
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
        profile.asEigen() += front.asEigen() * i->mapping;
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
    ndarray::Array<Pixel, 2, 2> mapping = ndarray::allocate(
        ndarray::makeVector(height, width)
    );
    ar >> boost::serialization::make_array(mapping.getData(), size);
    int constraintCount;
    ar >> constraintCount;
    CompoundShapeletBuilder builder(components, mapping);
    if (constraintCount > 0) {
        ndarray::Array<Pixel,2,2> cMatrix = ndarray::allocate(constraintCount, width);
        ndarray::Array<Pixel,1,1> cVector = ndarray::allocate(constraintCount);
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
    int height = _mapping.getSize<0>(), width = _mapping.getSize<1>();
    int size = width*height;
    ndarray::Array<Pixel, 2, 2> mapping = ndarray::copy(_mapping);
    ar << height;
    ar << width;
    ar << boost::serialization::make_array(mapping.getData(), size);
    int constraintCount = getConstraintCount();
    ar << constraintCount;
    if (constraintCount > 0) {
        ndarray::Array<Pixel,2,2> cMatrix = ndarray::copy(getConstraintMatrix());
        ndarray::Array<Pixel,1,1> cVector = ndarray::copy(getConstraintVector());
        ar << boost::serialization::make_array(cMatrix.getData(), cMatrix.getNumElements());
        ar << boost::serialization::make_array(cVector.getData(), cVector.getNumElements());
    }
}

//====================================== CompoundShapeletBuilder =========================================

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components,
    afwShapelets::BasisTypeEnum basisType,
    bool radialOnly
) :
    detail::CompoundShapeletBase(components)
{
    if (basisType == afwShapelets::HERMITE) {
        if (radialOnly) throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot specify radialOnly with HERMITE basis type."
        );
        return;
    }
    int offset = 0;
    std::vector<int> indices;
    for (ElementVector::iterator i = _elements.begin(); i != _elements.end(); ++i) {
        int size = i->component->getSize();
        ndarray::Array<Pixel,2,1> fBlock = _mapping[ndarray::view(offset, offset+size)(offset, offset+size)];
        for (int n = 0; n < i->component->getSize(); ++n) {
            afwShapelets::ConversionMatrix::convertCoefficientVector(
                fBlock[n], afwShapelets::HERMITE, afwShapelets::LAGUERRE, i->component->getOrder()
            );
        }
        if (radialOnly) {
            for (int n = 0; n <= i->component->getOrder(); n += 2) {
                indices.push_back(offset + afwShapelets::computeOffset(n) + n);
            }
        }
        offset += i->component->getSize();
    }
    if (radialOnly) {
        ndarray::Array<Pixel,2,2> newMapping(ndarray::allocate(_mapping.getSize<0>(), indices.size()));
        newMapping.deep() = 0.0;
        ndarray::Array<Pixel,2,0> newMappingT(newMapping.transpose());
        ndarray::Array<Pixel,2,0> oldMappingT(_mapping.transpose());
        for (int i = 0; i < newMappingT.getSize<0>(); ++i) {
            newMappingT[i].deep() = oldMappingT[indices[i]];
        }
        _mapping = newMapping;
        _resetElements();
    }
}

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components,
    ndarray::Array<Pixel const,2,1> const & mapping
) : 
    detail::CompoundShapeletBase(components, ndarray::copy(mapping)) 
{}

void CompoundShapeletBuilder::orthogonalize() {
    Eigen::MatrixXd v = computeInnerProductMatrix();
    Eigen::LLT<Eigen::MatrixXd> cholesky(v);
    Eigen::MatrixXd m = Eigen::MatrixXd::Identity(v.rows(), v.cols());
    cholesky.matrixL().transpose().solveInPlace(m);
    Matrix newMapping(ndarray::allocate(_mapping.getShape()));
    newMapping = _mapping.asEigen() * m;
    if (_constraintMatrix.getSize<0>() > 0) {
        ndarray::Array<Pixel,2,2> newConstraintMatrix = ndarray::allocate(_constraintMatrix.getShape());
        newConstraintMatrix.asEigen() = _constraintMatrix.asEigen() * m;
    }
    _mapping = newMapping.shallow();
    _resetElements();
}

void CompoundShapeletBuilder::normalizeFlux(int n) {
    ndarray::Array<Pixel,1,1> integral(ndarray::allocate(getSize()));
    integral.deep() = 0.0;
    integrate(integral);
    if (_constraintMatrix.getSize<0>() > 0) {
        _constraintMatrix = ndarray::copy(_constraintMatrix / integral[n]);
    }
    _mapping = ndarray::copy(_mapping / integral[n]);
    _resetElements();
}

void CompoundShapeletBuilder::slice(int start, int stop) {
    _mapping = _mapping[ndarray::view()(start, stop)];
    _constraintMatrix = _constraintMatrix[ndarray::view()(start, stop)];
    _resetElements();

}

void CompoundShapeletBuilder::setMapping(
    ndarray::Array<Pixel const,2,1> const & mapping
) {
    detail::checkSize(
        _mapping.getSize<0>(), mapping.getSize<0>(),
        "Aggregate component shapelet basis size (%d) does not match mapping mapping rows (%d)."
    );
    _mapping = ndarray::copy(mapping);
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

namespace {

class TruncatedDeVaucouleurProfileFunction : public ProfileFunction {
public:

    virtual double operator()(double radius) const {
        static double const DEVFAC = -7.66925;
        static double const DEVOUT = 8.0;
        static double const DEVCUT = 7.0;
        if (radius > DEVOUT) return 0.0;
        double z = 1.0;
        if (radius > DEVCUT) {
            z = (radius - DEVCUT) / (DEVOUT - DEVCUT);
            z = 1.0 - z*z;
        }
        z *= std::exp(DEVFAC * (std::pow(radius * radius + 0.0004, 0.125) - 1.0));
        return z;
    }

};

class TruncatedExponentialProfileFunction : public ProfileFunction {
public:

    double operator()(double radius) const {
        static double const EXPFAC = -1.67835;
        static double const EXPOUT = 4.0;
        static double const EXPCUT = 3.0;
        if (radius > EXPOUT) return 0.0;
        double z = 1.0;
        if (radius > EXPCUT) {
            z = (radius - EXPCUT) / (EXPOUT - EXPCUT);
            z = 1.0 - z*z;
        }
        z *= std::exp(EXPFAC * (radius - 1.0));
        return z;
    }

};

} // anonymous

ProfileFunction::Ptr ProfileFunction::makeTruncatedDeVaucouleur() {
    return ProfileFunction::Ptr(new TruncatedDeVaucouleurProfileFunction());
}

ProfileFunction::Ptr ProfileFunction::makeTruncatedExponential() {
    return ProfileFunction::Ptr(new TruncatedExponentialProfileFunction());
}

CompoundShapeletBuilder CompoundShapeletBuilder::approximate(
    ProfileFunction const & profile,
    ComponentVector const & components,
    double sersicRadius,
    double maxRadius,
    ndarray::Array<Pixel const,1,1> const & matchRadii
) {
    CompoundShapeletBuilder builder(components, afwShapelets::LAGUERRE, true);
    CompoundShapeletModelBasis::Ptr basis = builder.build();
    afw::geom::ellipses::Ellipse ellipse(EllipseCore(0.0, 0.0, 1.0));
    afw::detection::Footprint::Ptr footprint(new afw::detection::Footprint(afw::geom::Point2I(), maxRadius));
    ndarray::Array<Pixel,2,2> iConstraintMatrix(
        ndarray::allocate(footprint->getArea() + basis->getSize(), basis->getSize())
    );
    ndarray::Array<Pixel,2,2> eConstraintMatrix(
        ndarray::allocate(matchRadii.getSize<0>(), basis->getSize())
    );
    basis->evaluateRadialProfile(eConstraintMatrix, matchRadii);
    ndarray::Array<Pixel,1,1> eConstraintVector = ndarray::copy(matchRadii);
    for (
        ndarray::Array<Pixel,1,1>::Iterator i = eConstraintVector.begin();
        i != eConstraintVector.end();
        ++i
    ) {
        *i = profile(*i / sersicRadius);
    }
    iConstraintMatrix[ndarray::view(footprint->getArea(), iConstraintMatrix.getSize<0>())].asEigen().setIdentity();
    ndarray::Array<Pixel,2,2> modelMatrix(iConstraintMatrix[ndarray::view(0, footprint->getArea())]);
    basis->evaluate(modelMatrix, footprint, ellipse);
    ndarray::Array<Pixel,1,1> iConstraintVector(ndarray::allocate(iConstraintMatrix.getSize<0>()));
    ndarray::Array<Pixel,1,1> dataVector(ndarray::allocate(footprint->getArea()));
    iConstraintVector.deep() = 0.0;
    ndarray::Array<Pixel,1,1>::Iterator p = dataVector.begin();
    for (
        afw::detection::Footprint::SpanList::const_iterator i = footprint->getSpans().begin();
        i != footprint->getSpans().end();
        ++i
    ) {
        int y = (**i).getY();
        for (int x = (**i).getX0(); x <= (**i).getX1(); ++x) {
            double r = std::sqrt(x * x + y * y);
            *p = profile(r / sersicRadius);
            ++p;
        }
    }
    ndarray::Array<Pixel,1,1> rhs(ndarray::allocate(basis->getSize()));
    rhs.asEigen() = -(modelMatrix.asEigen().transpose() * dataVector.asEigen());
    ndarray::Array<Pixel,2,2> fisherMatrix(ndarray::allocate(basis->getSize(), basis->getSize()));
    fisherMatrix.asEigen()
        = modelMatrix.asEigen().transpose() * modelMatrix.asEigen();
    ndarray::Array<Pixel,1,1> coefficients(ndarray::allocate(basis->getSize()));

    double r = QPSolver(fisherMatrix, rhs)
        .inequality(iConstraintMatrix, iConstraintVector)
        .equality(eConstraintMatrix, eConstraintVector)
        .solve(coefficients);
    if (r == std::numeric_limits<double>::infinity()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "QP solver did not converge; program may be infeasible.");
    }

    ndarray::Array<Pixel,2,2> mapping(ndarray::allocate(builder.getMapping().getSize<0>(), 1));
    mapping.asEigen().col(0) = builder.getMapping().asEigen()
        * coefficients.asEigen();
    builder.setMapping(mapping);

    builder.dataImage = expandArray(*footprint, dataVector, footprint->getBBox());
    builder.modelImages = expandArray(*footprint, modelMatrix, footprint->getBBox());

    for (ElementVector::iterator i = builder._elements.begin(); i != builder._elements.end(); ++i) {
        i->component = ShapeletModelBasis::make(i->component->getOrder(), 
                                                i->component->getScale() / sersicRadius);
    }

    builder._constraintMatrix = ndarray::allocate(1, 1);
    builder._constraintMatrix[0] = 1.0;
    builder._constraintVector = ndarray::allocate(1);
    builder._constraintVector[0] = 0.0;

    return builder;
}



}}} //end namespace lsst::meas::multifit
