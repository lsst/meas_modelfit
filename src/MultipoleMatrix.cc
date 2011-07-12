#include "lsst/meas/multifit/MultipoleMatrix.h"
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

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

void MultipoleMatrix::applyMoments(
    lsst::afw::geom::ellipses::Ellipse & ellipse,
    lsst::ndarray::Array<Pixel const,1,1> const & coefficients
) const {
    Eigen::VectorXd m = ndarray::viewAsEigen(_array) * ndarray::viewAsEigen(coefficients);
    m.segment(1, 5) /= m[I0];
    m[IXX] -= m[IX] * m[IX];
    m[IYY] -= m[IY] * m[IY];
    m[IXY] -= m[IX] * m[IY];
    ellipse = afw::geom::ellipses::Ellipse(
        afw::geom::ellipses::Quadrupole(m[IXX], m[IYY], m[IXY]),
        afw::geom::Point2D(m[IX], m[IY])
    ).transform(ellipse.getGridTransform().invert());
}

}}} // namespace lsst::meas::multifit
