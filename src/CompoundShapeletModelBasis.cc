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
#include <Eigen/Cholesky>
#include <fstream>
#include "boost/serialization/binary_object.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

namespace afwShapelets = lsst::afw::math::shapelets;

namespace lsst { namespace meas { namespace multifit {

//========================================== ProfileFunction =========================================

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

//========================================== CompoundShapeletImpl =========================================

namespace detail {

// Ops that only work on ShapeletModelBasis instantiation of Impl go here.
class CompoundShapeletHelper {
public:

    static void convertBasis(
        CompoundShapeletImpl<ShapeletModelBasis> & impl,
        afwShapelets::BasisTypeEnum basisType,
        bool radialOnly
    );

    static Eigen::MatrixXd computeInnerProductMatrix(
        CompoundShapeletImpl<ShapeletModelBasis> const & impl
    );

    static void approximate(
        ProfileFunction const & profile,
        CompoundShapeletImpl<ShapeletModelBasis> & impl,
        double sersicRadius,
        double maxRadius,
        ndarray::Array<Pixel const,1,1> const & matchRadii
    );

};

template <typename Component>
class CompoundShapeletImpl {
public:
 
    typedef std::vector<PTR(Component)> ComponentVector;
    typedef typename std::vector<PTR(Component)>::const_iterator ComponentIter;
    typedef ndarray::EigenView<Pixel,2,1> Matrix;

    struct Element {
        PTR(Component) component;
        Matrix mapping;

        Element(
            PTR(Component) const & component_, 
            ndarray::Array<Pixel,2,1> const & fullMapping,
            int offset
        ) : 
            component(component_),
            mapping(fullMapping[ndarray::view(offset, offset + component->getSize())()])
        {}

        Element(PTR(Component) const & component_, Matrix const & mapping) :
            component(component_), mapping(mapping.getArray())
        {}

        Element & operator=(Element const & other) {
            if (&other != this) {
                component = other.component;
                mapping.setArray(other.mapping.getArray());
            }
            return *this;
        }
    };

    typedef std::vector<Element> ElementVector;
    typedef typename std::vector<Element>::const_iterator ElementIter;

    CompoundShapeletImpl(ComponentVector const & components, ndarray::Array<Pixel,2,1> const & mapping) :
        _elements(), _mapping(mapping)
    {
        checkSize(
            _computeSize(components), _mapping.getSize<0>(),
            "Aggregate component shapelet basis size (%d) does not match mapping mapping rows (%d)."
        );
        _fillElements(components);
    }

    CompoundShapeletImpl(
        ComponentVector const & components
    ) : _elements(), _mapping(_makeIdentity(_computeSize(components))) {
        _fillElements(components);
    }

    CompoundShapeletImpl(ElementVector & elements) : _elements(), _mapping() {
        _elements.swap(elements);
    }

    CompoundShapeletImpl(CompoundShapeletImpl const & other) :
        _elements(other._elements), _mapping(other._mapping)
    {}
    
    int getSize() const { return _mapping.getSize<1>(); }

    ComponentVector extractComponents() const {
        ComponentVector result;
        result.reserve(_elements.size());
        for (ElementIter i = _elements.begin(); i != _elements.end(); ++i) {
            result.push_back(i->component);
        }
        return result;
    }

    /// @brief Return the mapping from Shapelet components (rows) to compound basis functions (columns).
    lsst::ndarray::Array<Pixel const,2,1> getMapping() const { return _mapping; }

    void setMapping(lsst::ndarray::Array<Pixel,2,1> const & mapping) {
        _mapping = mapping;
        _resetElements();
    }

    ElementVector const & getElements() const { return _elements; }

    void slice(int start, int stop) {
        _mapping = _mapping[ndarray::view()(start, stop)];
        _resetElements();
    }

    void evaluate(
        lsst::ndarray::Array<Pixel,2,1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const {
        matrix.deep() = 0.0;
        for (ElementIter i = _elements.begin(); i != _elements.end(); ++i) {
            ndarray::Array<Pixel,2,2> front =
                ndarray::allocate(footprint->getArea(), i->component->getSize());
            i->component->evaluate(front, footprint, ellipse);
            ndarray::viewAsEigen(matrix) += ndarray::viewAsEigen(front) * i->mapping;
        }
    }

    void evaluateRadialProfile(
        ndarray::Array<Pixel,2,1> const & profile,
        ndarray::Array<Pixel const,1,1> const & radii
    ) const {
        profile.deep() = 0.0;
        for (ElementIter i = _elements.begin(); i != _elements.end(); ++i) {
            ndarray::Array<Pixel,2,2> front(ndarray::allocate(radii.getSize<0>(), i->component->getSize()));
            i->component->evaluateRadialProfile(front, radii);
            ndarray::viewAsEigen(profile) += ndarray::viewAsEigen(front) * i->mapping;
        }
    }

    void evaluateMultipoleMatrix(ndarray::Array<Pixel,2,2> const & matrix) const {
        matrix.deep() = 0.0;
        for (ElementIter i = _elements.begin(); i != _elements.end(); ++i) {
            ndarray::viewAsEigen(matrix) += 
                ndarray::viewAsEigen(i->component->getMultipoleMatrix().getArray()) * i->mapping;
        }
    }

private:

    friend class CompoundShapeletHelper;

    void _fillElements(ComponentVector const & components) {
        _elements.reserve(components.size());
        int offset = 0;
        for (ComponentIter i = components.begin(); i != components.end(); ++i) {
            _elements.push_back(Element(*i, _mapping, offset));
            offset += (**i).getSize();
        }
    }

    void _resetElements() {
        ElementVector new_elements;
        new_elements.reserve(_elements.size());
        int offset = 0; 
        for (ElementIter i = _elements.begin(); i != _elements.end(); ++i) {
            new_elements.push_back(Element(i->component, _mapping, offset));
            offset += i->component->getSize();
        }
        _elements.swap(new_elements);
    }

    static int _computeSize(ComponentVector const & components) {
        int size = 0;
        for (ComponentIter i = components.begin(); i != components.end(); ++i) {
            size += (**i).getSize();
        }
        return size;
    }

    static ndarray::Array<Pixel,2,2> _makeIdentity(int size) {
        ndarray::Array<Pixel,2,2> result(ndarray::allocate(size, size));
        ndarray::viewAsEigen(result).setIdentity();
        return result;
    }

    ElementVector _elements;
    ndarray::Array<Pixel,2,1> _mapping;
};

void CompoundShapeletHelper::convertBasis(
    CompoundShapeletImpl<ShapeletModelBasis> & impl,
    afwShapelets::BasisTypeEnum basisType,
    bool radialOnly
) {
    typedef CompoundShapeletImpl<ShapeletModelBasis> Impl;
    if (basisType == afwShapelets::HERMITE) {
        if (radialOnly) throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot specify radialOnly with HERMITE basis type."
        );
        return;
    }
    int offset = 0;
    std::vector<int> indices;
    for (Impl::ElementIter i = impl._elements.begin(); i != impl._elements.end(); ++i) {
        int size = i->component->getSize();
        ndarray::Array<Pixel,2,1> fBlock 
            = impl._mapping[ndarray::view(offset, offset+size)(offset, offset+size)];
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
        ndarray::Array<Pixel,2,2> newMapping(ndarray::allocate(impl._mapping.getSize<0>(), indices.size()));
        newMapping.deep() = 0.0;
        ndarray::Array<Pixel,2,0> newMappingT(newMapping.transpose());
        ndarray::Array<Pixel,2,0> oldMappingT(impl._mapping.transpose());
        for (int i = 0; i < newMappingT.getSize<0>(); ++i) {
            newMappingT[i].deep() = oldMappingT[indices[i]];
        }
        impl._mapping = newMapping;
        impl._resetElements();
    }
    
}

Eigen::MatrixXd CompoundShapeletHelper::computeInnerProductMatrix(
    CompoundShapeletImpl<ShapeletModelBasis> const & impl
) {
    typedef CompoundShapeletImpl<ShapeletModelBasis>::ElementVector::const_iterator Iter;
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(
        impl.getMapping().getSize<1>(), impl.getMapping().getSize<1>()
    );
    for (Iter i = impl.getElements().begin(); i != impl.getElements().end(); ++i) {
        for (Iter j = impl.getElements().begin(); j != impl.getElements().end(); ++j) {
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

void CompoundShapeletHelper::approximate(
    ProfileFunction const & profile,
    CompoundShapeletImpl<ShapeletModelBasis> & impl,
    double sersicRadius,
    double maxRadius,
    ndarray::Array<Pixel const,1,1> const & matchRadii
) {
    typedef CompoundShapeletImpl<ShapeletModelBasis> Impl;
    afw::geom::ellipses::Ellipse ellipse(EllipseCore(0.0, 0.0, 1.0));
    afw::detection::Footprint::Ptr footprint(new afw::detection::Footprint(afw::geom::Point2I(), maxRadius));
    ndarray::Array<Pixel,2,2> iConstraintMatrix(
        ndarray::allocate(footprint->getArea() + impl.getSize(), impl.getSize())
    );
    ndarray::Array<Pixel,2,2> eConstraintMatrix(
        ndarray::allocate(matchRadii.getSize<0>(), impl.getSize())
    );
    impl.evaluateRadialProfile(eConstraintMatrix, matchRadii);
    ndarray::Array<Pixel,1,1> eConstraintVector = ndarray::copy(matchRadii);
    for (
        ndarray::Array<Pixel,1,1>::Iterator i = eConstraintVector.begin();
        i != eConstraintVector.end();
        ++i
    ) {
        *i = profile(*i / sersicRadius);
    }
    ndarray::viewAsEigen(
        iConstraintMatrix[ndarray::view(footprint->getArea(), iConstraintMatrix.getSize<0>())]
    ).setIdentity();
    ndarray::Array<Pixel,2,2> modelMatrix(iConstraintMatrix[ndarray::view(0, footprint->getArea())]);
    impl.evaluate(modelMatrix, footprint, ellipse);
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
    ndarray::Array<Pixel,1,1> rhs(ndarray::allocate(impl.getSize()));
    ndarray::viewAsEigen(rhs) 
        = -(ndarray::viewAsTransposedEigen(modelMatrix) * ndarray::viewAsEigen(dataVector)).lazy();
    ndarray::Array<Pixel,2,2> fisherMatrix(ndarray::allocate(impl.getSize(), impl.getSize()));
    ndarray::viewAsEigen(fisherMatrix)
        = ndarray::viewAsTransposedEigen(modelMatrix) * ndarray::viewAsEigen(modelMatrix);
    ndarray::Array<Pixel,1,1> coefficients(ndarray::allocate(impl.getSize()));

    double r = QPSolver(fisherMatrix, rhs)
        .inequality(iConstraintMatrix, iConstraintVector)
        .equality(eConstraintMatrix, eConstraintVector)
        .solve(coefficients);
    if (r == std::numeric_limits<double>::infinity()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "QP solver did not converge; program may be infeasible.");
    }

    ndarray::Array<Pixel,2,2> mapping(ndarray::allocate(impl.getMapping().getSize<0>(), 1));
    ndarray::viewAsEigen(mapping).col(0) = ndarray::viewAsEigen(impl.getMapping())
        * ndarray::viewAsEigen(coefficients);
    impl.setMapping(mapping);

    ndarray::Array<Pixel,2,2> multipoleMatrix(ndarray::allocate(6, 1));
    impl.evaluateMultipoleMatrix(multipoleMatrix);
    double radius = std::sqrt(multipoleMatrix[MultipoleMatrix::IXX][0] 
                              / multipoleMatrix[MultipoleMatrix::I0][0]);
    for (Impl::ElementVector::iterator i = impl._elements.begin(); i != impl._elements.end(); ++i) {
        i->component = ShapeletModelBasis::make(
            i->component->getOrder(), i->component->getScale() / radius
        );
    }
}

} // namespace detail

//================================ ConvolvedCompoundShapeletModelBasis =====================================

namespace {

class ConvolvedCompoundShapeletModelBasis : public ModelBasis {
public:

    typedef detail::CompoundShapeletImpl<ModelBasis> Impl;

    ConvolvedCompoundShapeletModelBasis(
        int size, Impl::ElementVector & elements, ndarray::Array<Pixel const,2,2> const & multipoleMatrix
    ) : 
        ModelBasis(size), _impl(new Impl(elements))
    {
        attachMultipoleMatrix(multipoleMatrix);
    }

protected:

    virtual void _evaluate(
        lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        lsst::afw::geom::Ellipse const & ellipse
    ) const {
        _impl->evaluate(matrix, footprint, ellipse);
    }

private:
    boost::scoped_ptr<Impl> _impl;
};

} // anonymous

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
    typedef ConvolvedCompoundShapeletModelBasis::Impl ConvolvedImpl;
    ConvolvedImpl::ElementVector convolvedElements;
    convolvedElements.reserve(_impl->getElements().size());
    for (Impl::ElementIter i = _impl->getElements().begin(); i != _impl->getElements().end(); ++i) {
        convolvedElements.push_back(
            ConvolvedImpl::Element(i->component->convolve(psf), i->mapping)
        );
    }
    return ModelBasis::Ptr(
        new ConvolvedCompoundShapeletModelBasis(
            this->getSize(), boost::ref(convolvedElements), getMultipoleMatrix().getArray()
        )
    );
}

ModelBasis::Ptr CompoundShapeletModelBasis::convolve(
    afwShapelets::MultiShapeletFunction const & psf
) const {
    typedef ConvolvedCompoundShapeletModelBasis::Impl ConvolvedImpl;
    ConvolvedImpl::ElementVector convolvedElements;
    convolvedElements.reserve(_impl->getElements().size() * psf.getElements().size());
    for (Impl::ElementIter i = _impl->getElements().begin(); i != _impl->getElements().end(); ++i) {
        for (LocalPsf::MultiShapelet::ElementList::const_iterator j = psf.getElements().begin(); 
             j != psf.getElements().end();
             ++j
        ) {
            convolvedElements.push_back(
                ConvolvedImpl::Element(i->component->convolve(*j), i->mapping)
            );
        }
    }
    return ModelBasis::Ptr(
        new ConvolvedCompoundShapeletModelBasis(
            this->getSize(), boost::ref(convolvedElements), getMultipoleMatrix().getArray()
        )
    );
}
    
CompoundShapeletModelBasis::ComponentVector CompoundShapeletModelBasis::extractComponents() const {
    return _impl->extractComponents();
}

lsst::ndarray::Array<Pixel const,2,1> CompoundShapeletModelBasis::getMapping() const {
    return _impl->getMapping();
}

Eigen::MatrixXd CompoundShapeletModelBasis::computeInnerProductMatrix() const {
    return detail::CompoundShapeletHelper::computeInnerProductMatrix(*_impl);
}

CompoundShapeletModelBasis::~CompoundShapeletModelBasis() {}

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
    for (int i =0; i < nComponents; ++i) {
        ar >> order;
        ar >> scale;
        components[i] = multifit::ShapeletModelBasis::make(order, scale);            
    }
    int width, height;
    ar >> height;
    ar >> width;
    int size = width*height;
    ndarray::Array<Pixel,2,2> mapping = ndarray::allocate(
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

    int nElement = _impl->getElements().size();
    ar << nElement;
    for (Impl::ElementIter i = _impl->getElements().begin(); i != _impl->getElements().end(); ++i) {
        int order = i->component->getOrder();
        double scale = i->component->getScale();
        ar << order;
        ar << scale;
    }
    ndarray::Array<Pixel,2,2> mapping = ndarray::copy(_impl->getMapping());
    int height = mapping.getSize<0>(), width = mapping.getSize<1>();
    int size = width*height;
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

void CompoundShapeletModelBasis::_evaluate(
    ndarray::Array<Pixel, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
    afw::geom::Ellipse const & ellipse
) const {
    _impl->evaluate(matrix, footprint, ellipse);
}

void CompoundShapeletModelBasis::_evaluateRadialProfile(
    lsst::ndarray::Array<Pixel,2,1> const & profile,
    lsst::ndarray::Array<Pixel const,1,1> const & radii
) const {
    _impl->evaluateRadialProfile(profile, radii);
}

CompoundShapeletModelBasis::CompoundShapeletModelBasis(
    boost::shared_ptr<Impl> const & impl,
    ndarray::Array<Pixel const,2,2> const & multipoleMatrix,
    ndarray::Array<Pixel const,2,1> const & constraintMatrix,
    ndarray::Array<Pixel const,1,1> const & constraintVector
) : ModelBasis(impl->getSize()), _impl(impl) {
    if (constraintMatrix.getSize<0>() > 0) {
        attachConstraint(constraintMatrix, constraintVector);
    }
    attachMultipoleMatrix(multipoleMatrix);
}

//====================================== CompoundShapeletBuilder =========================================

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components,
    afwShapelets::BasisTypeEnum basisType,
    bool radialOnly
) : _impl(new Impl(components)) {}

CompoundShapeletBuilder::CompoundShapeletBuilder(
    ComponentVector const & components,
    ndarray::Array<Pixel const,2,1> const & mapping
) : _impl(new Impl(components, ndarray::copy(mapping))) {}

CompoundShapeletBuilder CompoundShapeletBuilder::approximate(
    ProfileFunction const & profile,
    ComponentVector const & components,
    double sersicRadius,
    double maxRadius,
    ndarray::Array<Pixel const,1,1> const & matchRadii
) {
    boost::shared_ptr<Impl> impl(new Impl(components));
    detail::CompoundShapeletHelper::convertBasis(*impl, afwShapelets::LAGUERRE, true);
    detail::CompoundShapeletHelper::approximate(profile, *impl, sersicRadius, maxRadius, matchRadii);
    CompoundShapeletBuilder builder(impl);
    builder._constraintMatrix = ndarray::allocate(1, 1);
    builder._constraintMatrix[0] = 1.0;
    builder._constraintVector = ndarray::allocate(1);
    builder._constraintVector[0] = 0.0;
    return builder;
}

int CompoundShapeletBuilder::getSize() const { return _impl->getSize(); }

void CompoundShapeletBuilder::normalizeFlux(int n) {
    if (!_impl.unique()) {
        boost::shared_ptr<Impl> newImpl(new Impl(*_impl));
        newImpl.swap(_impl);
    }
    ndarray::Array<Pixel,1,1> integration(ndarray::allocate(getSize()));
    evaluateIntegration(integration);
    if (_constraintMatrix.getSize<0>() > 0) {
        _constraintMatrix = ndarray::copy(_constraintMatrix / integration[n]);
    }
    ndarray::Array<Pixel,2,1> newMapping(ndarray::copy(_impl->getMapping()));
    newMapping.deep() /= integration[n];
    _impl->setMapping(newMapping);
}

void CompoundShapeletBuilder::orthogonalize() {
    if (!_impl.unique()) {
        boost::shared_ptr<Impl> newImpl(new Impl(*_impl));
        newImpl.swap(_impl);
    }
    Eigen::MatrixXd v = detail::CompoundShapeletHelper::computeInnerProductMatrix(*_impl);
    Eigen::LLT<Eigen::MatrixXd> cholesky(v);
    Eigen::MatrixXd m = Eigen::MatrixXd::Identity(v.rows(), v.cols());
    cholesky.matrixL().transpose().solveTriangularInPlace(m);
    ndarray::Array<Pixel,2,2> newMapping(ndarray::allocate(_impl->getMapping().getShape()));
    ndarray::viewAsEigen(newMapping) = ndarray::viewAsEigen(_impl->getMapping()) * m;
    if (_constraintMatrix.getSize<0>() > 0) {
        ndarray::Array<Pixel,2,2> newConstraintMatrix = ndarray::allocate(_constraintMatrix.getShape());
        ndarray::viewAsEigen(newConstraintMatrix) = ndarray::viewAsEigen(_constraintMatrix) * m;
    }
    _impl->setMapping(newMapping);
}

void CompoundShapeletBuilder::slice(int start, int stop) {
    if (!_impl.unique()) {
        boost::shared_ptr<Impl> newImpl(new Impl(*_impl));
        newImpl.swap(_impl);
    }
    _constraintMatrix = _constraintMatrix[ndarray::view()(start, stop)];
    _impl->slice(start, stop);
}

CompoundShapeletBuilder::ComponentVector CompoundShapeletBuilder::extractComponents() const {
    return _impl->extractComponents();
}

ndarray::Array<Pixel const,2,1> CompoundShapeletBuilder::getMapping() const {
    return _impl->getMapping();
}

void CompoundShapeletBuilder::setMapping(
    ndarray::Array<Pixel const,2,1> const & mapping
) {
    detail::checkSize(
        _impl->getMapping().getSize<0>(), mapping.getSize<0>(),
        "Aggregate component shapelet basis size (%d) does not match mapping mapping rows (%d)."
    );
    if (!_impl.unique()) {
        boost::shared_ptr<Impl> newImpl(new Impl(*_impl));
        newImpl.swap(_impl);
    }
    _impl->setMapping(ndarray::copy(mapping));
}

Eigen::MatrixXd CompoundShapeletBuilder::computeInnerProductMatrix() const {
    return detail::CompoundShapeletHelper::computeInnerProductMatrix(*_impl);
}

void CompoundShapeletBuilder::evaluate(
    lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
    lsst::afw::geom::Ellipse const & ellipse
) const {
    _impl->evaluate(matrix, footprint, ellipse);
}

void CompoundShapeletBuilder::evaluateRadialProfile(
    lsst::ndarray::Array<Pixel,2,1> const & profile,
    lsst::ndarray::Array<Pixel const,1,1> const & radii
) const {
    _impl->evaluateRadialProfile(profile, radii);
}

void CompoundShapeletBuilder::evaluateIntegration(lsst::ndarray::Array<Pixel,1,1> const & vector) const {
    vector.deep() = 0.0;
    for (Impl::ElementIter i = _impl->getElements().begin(); i != _impl->getElements().end(); ++i) {
        ndarray::viewAsTransposedEigen(vector) += 
            ndarray::viewAsTransposedEigen(i->component->getIntegration()) * i->mapping;
    }
}

void CompoundShapeletBuilder::evaluateMultipoleMatrix(lsst::ndarray::Array<Pixel,2,2> const & matrix) const {
    _impl->evaluateMultipoleMatrix(matrix);
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
    ndarray::Array<Pixel,2,2> multipoleMatrix(ndarray::allocate(6, _impl->getSize()));
    evaluateMultipoleMatrix(multipoleMatrix);
    return CompoundShapeletModelBasis::Ptr(
        new CompoundShapeletModelBasis(_impl, multipoleMatrix, _constraintMatrix, _constraintVector)
    );
}

CompoundShapeletBuilder::~CompoundShapeletBuilder() {}

}}} //end namespace lsst::meas::multifit
