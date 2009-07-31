#ifndef LSST_MEAS_MULTIFIT_AFFINETRANSFORM_H
#define LSST_MEAS_MULTIFIT_AFFINETRANSFORM_H

#include "Eigen/Core"

namespace Eigen {
    typedef Vector<double, 6, 1> Vector6d;
}

namespace lsst {
namespace meas {
namespace multifit {        


typedef Eigen::Transform2d TransformType;

class AffineTransform {
    enum VectorIndices {XX=0, YX=1, XY=2, YY=3, X=4, Y=5};

    AffineTransform() : _transform(Eigen::Matrix3d::Identity()) {}
    explicit AffineTransform(AffineTransform const & other) 
        : _transform(other._transform) {}
    explicit AffineTransform(Eigen::Matrix3d const & matrix)
        : _transform(matrix) {}
    explicit AffineTransform(Eigen::Matrix2d const & m, Eigen::Vector2d const v)
        : _transform(Eigen::Translation2d(v)*TransformType(m)) {}
    explicit AffineTransform(TransformType const & transform) 
        : _transform(transform) {}
    explicit AffineTransform(Vector6d const & vector) {
        set(vector);
    }
    AffineTransform invert() {
        Eigen::LU<Eigen::Matrix2d> lu(_transform.linear());
        if (!lu.isInvertible()) throw SingularTransformError();
        Eigen::Matrix2d inv = lu.inverse();
        return AffineTransform(inv,-inv*_transform.translation());  
    }
    static AffineTransform makeScaling(double const s) {
        return AffineTransform(TransformType(Eigen::Scaling<double, 2>(s)));
    }
    void set(TransformType const & t) {_transform = t}
    void set(Eigen::Matrix3d const & m) {
        _transform = m;    
    }
    void set(Vector6d const & v) {
        _transform.matrix().col(0) << v.start<2>(), 0;
        _transform.matrix().col(1) << v.segment<2>(2), 0;
        _transform.matrix().col(2) << v.end<2>(), 1;
    }
    void fillVector(Vector6d & v) const {
        v << _transform.matrix().col(0).start<2>(), 
                _transform.matrix().col(1).start<2>(),
                _transform.matrix().col(2).start<2>();    
    }
    Vector6d getVector() {
        return ((Vector6d() <<_transform.matrix().col(0).start<2>(), 
                _transform.matrix().col(1).start<2>(),
                _transform.matrix().col(2).start<2>()).finished);    
    }

    TransformType & getTransform() {return _transform;}

    friend AffineTransform operator *(AffineTransform const & a,
            AffineTransform const & b) {
        return AffineTransform(a._transform * b._transform);    
    }
    double operator[](int const & i) {
        _return transform.data()[_MATRIX_INDEX[i]];
    }
protected:
    TransformType _transform;    

    static int[6] _MATRIX_INDEX = {0, 1, 3, 4, 6, 7};
}

}}} //namespace lsst::meas::multifit
#endif //LSST_MEAS_MULTIFIT_AFFINETRANSFORM_H

