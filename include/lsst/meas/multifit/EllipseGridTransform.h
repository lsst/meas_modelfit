#ifndef LSST_MEAS_MULTIFIT_ELLIPSE_GRID_TRANSFORM_H
#define LSST_MEAS_MULTIFIT_ELLIPSE_GRID_TRANSFORM_H

#include <Eigen/Core>
#include "lsst/afw/geom.h"


namespace lsst {
namespace meas {
namespace multifit {

/**
 *  A differentiable expression object that computes the LinearTransform that, 
 *  when applied to a Fourier-space coordinate grid, maps the unit circle onto 
 *  an ellipse in real-space.
 */
class EllipseGridTransform {
public:
    typedef boost::shared_ptr<EllipseGridTransform> Ptr;
    typedef boost::shared_ptr<EllipseGridTransform const> ConstPtr;

    typedef Eigen::Matrix<double,4,3> DerivativeMatrix;
    
    explicit EllipseGridTransform(
        lsst::afw::geom::ellipses::BaseCore const & ellipse,
        lsst::afw::geom::ExtentI const & realDimensions
    );

    operator lsst::afw::geom::LinearTransform () const {
        return lsst::afw::geom::LinearTransform(
            _matrices->scaling * _matrices->tail
        );
    }

    DerivativeMatrix dEllipse() const;

private:
    struct Matrices {
        Eigen::Matrix2d factor;
        Eigen::Matrix2d scaling;
        Eigen::Matrix2d rotation;
        Eigen::Matrix2d tail;
        Eigen::Matrix3d jacobian;
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW; 
    };

    lsst::afw::geom::ellipses::Axes _axes;
    boost::shared_ptr<Matrices const> _matrices;

    static Eigen::Matrix2d const & _dA();
    static Eigen::Matrix2d const & _dB();
};

}}} // namespace lsst::meas::multifit

#endif 
