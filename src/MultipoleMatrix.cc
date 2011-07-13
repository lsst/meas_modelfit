#include "lsst/meas/multifit/MultipoleMatrix.h"
#include "lsst/meas/multifit/grid/Grid.h"
#include "lsst/ndarray/eigen.h"


namespace lsst { namespace meas { namespace multifit {

#if 0
void MultipoleMatrix::transform(afw::geom::AffineTransform const & t) {
    typedef afw::geom::AffineTransform AT;
    ndarray::Array<Pixel,2,2> old(ndarray::copy(_array));
    _array[IX] = t[AT::XX] * old[IX] + t[AT::XY] * old[IY];
    _array[IY] = t[AT::YX] * old[IX] + t[AT::YY] * old[IY];
    _array[IXX] = t[AT::XX] * t[AT::XX] * old[IXX] + t[AT::XY] * t[AT::XY] * old[IYY]
        + 2.0 * t[AT::XX] * t[AT::XY] * old[IXY];
    _array[IYY] = t[AT::YX] * t[AT::YX] * old[IXX] + t[AT::YY] * t[AT::YY] * old[IYY]
        + 2.0 * t[AT::YX] * t[AT::YY] * old[IXY];
    _array[IXY] = t[AT::XX] * t[AT::YX] * old[IXX] + t[AT::XY] * t[AT::YY] * old[IYY]
        + (t[AT::XY] * t[AT::YX] + t[AT::XX] * t[AT::YY]) * old[IXY];
}

void MultipoleMatrix::scale(double factor) {
    _array[IX] *= factor;
    _array[IY] *= factor;
    _array[IXX] *= factor * factor;
    _array[IYY] *= factor * factor;
    _array[IXY] *= factor * factor;
}
#endif

lsst::afw::geom::Point2D MultipoleMatrix::computeCentroid(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) const {
    return ellipse.getGridTransform().invert()(computeCentroid(coefficients));
}

lsst::afw::geom::Point2D MultipoleMatrix::computeCentroid(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients
) const {
    double m0 = ndarray::viewAsEigen(_array[I0]).dot(ndarray::viewAsEigen(coefficients));
    double mx = ndarray::viewAsEigen(_array[IX]).dot(ndarray::viewAsEigen(coefficients));
    double my = ndarray::viewAsEigen(_array[IY]).dot(ndarray::viewAsEigen(coefficients));
    return afw::geom::Point2D(mx / m0, my / m0);
}

lsst::afw::geom::ellipses::Quadrupole MultipoleMatrix::computeQuadrupole(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
    lsst::afw::geom::ellipses::Ellipse const & ellipse,
    lsst::afw::geom::Point2D const & centroid
) const {
    afw::geom::AffineTransform transform = ellipse.getGridTransform();
    return computeQuadrupole(coefficients, transform(centroid)).transform(transform.getLinear().invert());
}

lsst::afw::geom::ellipses::Quadrupole MultipoleMatrix::computeQuadrupole(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
    lsst::afw::geom::Point2D const & centroid
) const {
    double m0 = ndarray::viewAsEigen(_array[I0]).dot(ndarray::viewAsEigen(coefficients));
    double mxx = ndarray::viewAsEigen(_array[IXX]).dot(ndarray::viewAsEigen(coefficients));
    double myy = ndarray::viewAsEigen(_array[IYY]).dot(ndarray::viewAsEigen(coefficients));
    double mxy = ndarray::viewAsEigen(_array[IXY]).dot(ndarray::viewAsEigen(coefficients));
    mxx /= m0;
    myy /= m0;
    mxy /= m0;
    mxx -= centroid.getX() * centroid.getX();
    myy -= centroid.getY() * centroid.getY();
    mxy -= centroid.getX() * centroid.getY();
    return afw::geom::ellipses::Quadrupole(mxx, myy, mxy);
}

lsst::afw::geom::ellipses::Ellipse MultipoleMatrix::computeEllipse(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients,
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) const {
    return computeEllipse(coefficients).transform(ellipse.getGridTransform().invert());
}

lsst::afw::geom::ellipses::Ellipse MultipoleMatrix::computeEllipse(
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients
) const {
    Eigen::VectorXd m = ndarray::viewAsEigen(_array) * ndarray::viewAsEigen(coefficients);
    m.segment(1, 5) /= m[I0];
    m[IXX] -= m[IX] * m[IX];
    m[IYY] -= m[IY] * m[IY];
    m[IXY] -= m[IX] * m[IY];
    return afw::geom::ellipses::Ellipse(
        afw::geom::ellipses::Quadrupole(m[IXX], m[IYY], m[IXY]),
        afw::geom::Point2D(m[IX], m[IY])
    );
}

void MultipoleMatrix::updateParameters(
    Grid::Ptr const & grid, 
    lsst::ndarray::Array<double const,1,1> const & parametersIn,
    lsst::ndarray::Array<double,1,1> const & parametersOut,
    lsst::ndarray::Array<Pixel const,1,1> & coefficients
) {
    typedef afw::geom::AffineTransform AT;
    ndarray::Array<Pixel,2,2> pointSourceMultipoleMatrixArray(ndarray::allocate(6, 1));
    pointSourceMultipoleMatrixArray.deep() = 0.0;
    pointSourceMultipoleMatrixArray[0] = 1.0;
    ndarray::EigenView<Pixel const,1,1> c(coefficients);
    parametersOut.deep() = parametersIn;
    for (
        Grid::PositionArray::iterator positionIter = grid->positions.begin();
        positionIter != grid->positions.end();
        ++positionIter
    ) {
        double m0 = 0.0, mx = 0.0, my = 0.0;
        for (
            Grid::ObjectComponentArray::iterator objectIter = grid->objects.begin();
            objectIter != grid->objects.end();
            ++objectIter
        ) {
            if (objectIter->getPosition() != grid::PositionElement::Ptr(positionIter)) continue;
            ndarray::EigenView<Pixel const,2,2> mm(pointSourceMultipoleMatrixArray);
            afw::geom::AffineTransform t;
            afw::geom::ellipses::Ellipse ellipse = objectIter->makeEllipse(parametersIn);
            if (objectIter->getBasis()) {
                mm.setArray(objectIter->getBasis()->getMultipoleMatrix().getArray());
                t = ellipse.getGridTransform().invert();
            }
            double s0 = 0.0, sx = 0.0, sy = 0.0;
            int n = objectIter->getFluxGroup()->isVariable() ? grid->frames.size() : grid->getFilterCount();
            for (int i = 0; i < n; ++i) {
                int offset = objectIter->getFluxGroup()->getCoefficientOffset(i)
                    + objectIter->getGroupCoefficientOffset();
                s0 += mm.row(I0).dot(c.segment(offset, objectIter->getSourceCoefficientCount()));
                sx += mm.row(IX).dot(c.segment(offset, objectIter->getSourceCoefficientCount()));
                sy += mm.row(IY).dot(c.segment(offset, objectIter->getSourceCoefficientCount()));
            }
            m0 += s0;
            mx += t[AT::XX] * sx + t[AT::XY] * sy + t[AT::X] * s0;
            my += t[AT::YX] * sx + t[AT::YY] * sy + t[AT::Y] * s0;
        }
        parametersOut[positionIter->offset + 0] += mx / m0;
        parametersOut[positionIter->offset + 1] += my / m0;
    }
    for (
        Grid::RadiusArray::iterator radiusIter = grid->radii.begin();
        radiusIter != grid->radii.end();
        ++radiusIter
    ) {
        double m0 = 0.0, mxx = 0.0, myy = 0.0, mxy = 0.0;
        for (
            Grid::ObjectComponentArray::iterator objectIter = grid->objects.begin();
            objectIter != grid->objects.end();
            ++objectIter
        ) {
            if (objectIter->getRadius() != grid::RadiusElement::Ptr(radiusIter)) continue;
            ndarray::EigenView<Pixel const,2,2> mm(objectIter->getBasis()->getMultipoleMatrix().getArray());
            afw::geom::ellipses::Ellipse ellipse = objectIter->makeEllipse(parametersIn);
            afw::geom::Point2D p = objectIter->makePoint(parametersOut);
            afw::geom::AffineTransform t = ellipse.getGridTransform().invert();
            Eigen::VectorXd s = Eigen::VectorXd::Zero(6);
            int n = objectIter->getFluxGroup()->isVariable() ? grid->frames.size() : grid->getFilterCount();
            for (int i = 0; i < n; ++i) {
                int offset = objectIter->getFluxGroup()->getCoefficientOffset(i)
                    + objectIter->getGroupCoefficientOffset();
                s += mm * c.segment(offset, objectIter->getSourceCoefficientCount());
            }
            m0 += s[I0];
            mxx += t[AT::XX] * t[AT::XX] * s[IXX] + t[AT::XY] * t[AT::XY] * s[IYY] 
                + 2.0 * t[AT::XX] * t[AT::XY] * s[IXY]
                + 2.0 * t[AT::XX] * (t[AT::X] - p.getX()) * s[IX]
                + 2.0 * t[AT::XY] * (t[AT::X] - p.getX()) * s[IY]
                + (t[AT::X] - p.getX()) * (t[AT::X] - p.getX()) * s[I0];
            myy += t[AT::YX] * t[AT::YX] * s[IXX] + t[AT::YY] * t[AT::YY] * s[IYY]
                + 2.0 * t[AT::YX] * t[AT::YY] * s[IXY]
                + 2.0 * t[AT::YX] * (t[AT::Y] - p.getY()) * s[IX]
                + 2.0 * t[AT::YY] * (t[AT::Y] - p.getY()) * s[IY]
                + (t[AT::Y] - p.getY()) * (t[AT::Y] - p.getY()) * s[I0];
            mxy += t[AT::XX] * t[AT::YX] * s[IXX] + t[AT::XY] * t[AT::YY] * s[IYY]
                + (t[AT::XX] * t[AT::YY] + t[AT::XY] * t[AT::YX]) * s[IXY]
                + (t[AT::XX] * (t[AT::Y] - p.getY()) + t[AT::YX] * (t[AT::X] - p.getX())) * s[IX]
                + (t[AT::XY] * (t[AT::Y] - p.getY()) + t[AT::YY] * (t[AT::X] - p.getX())) * s[IY]
                + (t[AT::X] - p.getX()) * (t[AT::Y] - p.getY()) * s[I0];
        }
        EllipseCore core(afw::geom::ellipses::Quadrupole(mxx / m0, myy / m0, mxy / m0));
        parametersOut[radiusIter->offset] = std::sqrt(0.5*((mxx + myy) / m0));
    }
    for (
        Grid::EllipticityArray::iterator ellipticityIter = grid->ellipticities.begin();
        ellipticityIter != grid->ellipticities.end();
        ++ellipticityIter
    ) {
        double m0 = 0.0, mxx = 0.0, myy = 0.0, mxy = 0.0;
        for (
            Grid::ObjectComponentArray::iterator objectIter = grid->objects.begin();
            objectIter != grid->objects.end();
            ++objectIter
        ) {
            if (objectIter->getEllipticity() != grid::EllipticityElement::Ptr(ellipticityIter)) continue;
            ndarray::EigenView<Pixel const,2,2> mm(objectIter->getBasis()->getMultipoleMatrix().getArray());
            afw::geom::ellipses::Ellipse ellipse = objectIter->makeEllipse(parametersIn);
            afw::geom::Point2D p = objectIter->makePoint(parametersOut);
            afw::geom::AffineTransform t = ellipse.getGridTransform().invert();
            Eigen::VectorXd s = Eigen::VectorXd::Zero(6);
            int n = objectIter->getFluxGroup()->isVariable() ? grid->frames.size() : grid->getFilterCount();
            for (int i = 0; i < n; ++i) {
                int offset = objectIter->getFluxGroup()->getCoefficientOffset(i)
                    + objectIter->getGroupCoefficientOffset();
                s += mm * c.segment(offset, objectIter->getSourceCoefficientCount());
            }
            m0 += s[I0];
            mxx += t[AT::XX] * t[AT::XX] * s[IXX] + t[AT::XY] * t[AT::XY] * s[IYY] 
                + 2.0 * t[AT::XX] * t[AT::XY] * s[IXY]
                + 2.0 * t[AT::XX] * (t[AT::X] - p.getX()) * s[IX]
                + 2.0 * t[AT::XY] * (t[AT::X] - p.getX()) * s[IY]
                + (t[AT::X] - p.getX()) * (t[AT::X] - p.getX()) * s[I0];
            myy += t[AT::YX] * t[AT::YX] * s[IXX] + t[AT::YY] * t[AT::YY] * s[IYY]
                + 2.0 * t[AT::YX] * t[AT::YY] * s[IXY]
                + 2.0 * t[AT::YX] * (t[AT::Y] - p.getY()) * s[IX]
                + 2.0 * t[AT::YY] * (t[AT::Y] - p.getY()) * s[IY]
                + (t[AT::Y] - p.getY()) * (t[AT::Y] - p.getY()) * s[I0];
            mxy += t[AT::XX] * t[AT::YX] * s[IXX] + t[AT::XY] * t[AT::YY] * s[IYY]
                + (t[AT::XX] * t[AT::YY] + t[AT::XY] * t[AT::YX]) * s[IXY]
                + (t[AT::XX] * (t[AT::Y] - p.getY()) + t[AT::YX] * (t[AT::X] - p.getX())) * s[IX]
                + (t[AT::XY] * (t[AT::Y] - p.getY()) + t[AT::YY] * (t[AT::X] - p.getX())) * s[IY]
                + (t[AT::X] - p.getX()) * (t[AT::Y] - p.getY()) * s[I0];
        }
        EllipseCore core(afw::geom::ellipses::Quadrupole(mxx / m0, myy / m0, mxy / m0));
        parametersOut[ellipticityIter->offset + 0] = core.getE1();
        parametersOut[ellipticityIter->offset + 1] = core.getE2();
    }
}

}}} // namespace lsst::meas::multifit
