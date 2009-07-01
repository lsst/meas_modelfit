#ifndef LSST_MEAS_MULTIFIT_ELLIPSE_H
#define LSST_MEAS_MULTIFIT_ELLIPSE_H

#include "Eigen/Core"

namespace lsst {
namespace meas {
namespace multifit {

class Ellipse {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    typedef Eigen::Vector2d Coordinate;
    typedef Eigen::Matrix2d MomentsMatrix;
   
    class Axes;
    class Ellipticity;

    class Moments {
        Eigen::Vector3d _data;
    public:

        void setXX(double xx) { _data.x() = xx; }
        double getXX() const { return _data.x(); }

        void setYY(double yy) { _data.y() = yy; }
        double getYY() const { return _data.y(); }

        void setXY(double xy) { _data.z() = xy; }
        double getXY() const { return _data.z(); }

        double getDeterminant() const {
            return getXX()*getYY() - getXY()*getXY(); 
        }

        MomentsMatrix getMatrix() const {
            MomentsMatrix r; r << _data.x(), _data.z(), _data.y(), _data.z(); 
            return r;
        }

        bool isValid() const { return getDeterminant() >= 0; }

        template <typename Derived>
        explicit Moments(Eigen::MatrixBase<Derived> const & data) 
            : _data(data) {
        }

        explicit Moments(double xx=0, double yy=0, double xy=0)
            : _data(xx,yy,xy) {}

        Moments(Axes const & axes);
        Moments(Ellipticity const & ellipticity);
        
        template <typename Derived>
        void fill(Eigen::MatrixBase<Derived> & vector) const { vector = _data; }

        template <typename Derived>
        Moments& operator=(Eigen::MatrixBase<Derived> const & vector) { 
            _data = vector; 
            return *this; 
        }

        Moments& operator=(Axes const & axes);
        Moments& operator=(Ellipticity const & ellipticity);

        /** 
         *  Functor for use in evaluating elliptically symmetric Functions
         *  
         *
         *  RadialFraction computes, for a given point p, the squared ratio z
         *  of the radius of the ellipse (with zero centroid) along the ray 
         *  from the origin to p, divided by the distance from the origin to p.
         *  Evaluating a radial profile f(r) with unit radius with r = sqrt(z)
         *  will produce an elliptically symmetric function with that radial 
         *  profile and centroid, size, and ellipticity matching the ellipse.
         *
         *  For an ellipse with non-zero center, simply use let 
         *  p = p - ellipse.getCenter().
         */
        class RadialFraction {
        public:

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            RadialFraction(Moments const & moments) 
                : _matrix(moments.getMatrix()),
                  _inv_matrix(_matrix.inverse()) {}
	
            /** \brief Evaluate the RadialFraction z at the given point. */
            double operator()(Coordinate const & p) const {
                return std::sqrt(p.dot(_inv_matrix * p));
            }
            
            /** 
             *  \brief Evaluate the gradient of RadialFraction (derivative with
             *  respect to p).
             *
             *  The derivative with respect to the center of the ellipse is the
             *  negative gradient.
             */
            Eigen::RowVector2d computeGradient(Coordinate const & p) const {
                return Eigen::RowVector2d(-2.0*_inv_matrix*p);
            }

            /** \brief Evaluate the Jacobian of RadialFraction (derivative with
             *  respect to ellipse moments {xx,yy,xy}) at the given point. 
             */
            Eigen::RowVector3d computeJacobian(Coordinate const & p) const {
                Eigen::RowVector3d vec;
                Coordinate tmp1 = _matrix * p;
                MomentsMatrix tmp2;                
                tmp2.part<Eigen::SelfAdjoint>() = 
                        (tmp1 * tmp1.adjoint()).lazy();
                vec[0] = -tmp2(0,0);
                vec[1] = -tmp2(1,1);
                vec[2] = -2.0*tmp2(0,1);
                return vec;
            }

        protected:
	        MomentsMatrix _matrix;
            MomentsMatrix _inv_matrix;
        };
    };

    class Ellipticity {
        Eigen::Vector3d _data;
    public:

        void setE1(double e1) { _data.x() = e1; }
        double getE1() const { return _data.x(); }

        void setE2(double e2) { _data.y() = e2; }
        double getE2() const { return _data.y(); }

        void setRadius(double r) { _data.z() = r; }
        double getRadius() const { return _data.z(); }

        bool isValid() const { 
            return (getE1()*getE1() + getE2()*getE2()) <= 1 && getRadius() >= 0; 
        }

        template <typename Derived>
        Ellipticity(Eigen::MatrixBase<Derived> const & data) : _data(data) {}

        explicit Ellipticity(double e1=0, double e2=0, double radius=0) 
                : _data(e1,e2,radius) {
        }

        Ellipticity(Moments const & moments);
        Ellipticity(Axes const & axes);

        template <typename Derived>
        void fill(Eigen::MatrixBase<Derived> & vector) const { vector = _data; }

        template <typename Derived>
        Ellipticity& operator=(Eigen::MatrixBase<Derived> const & vec) { 
            _data = vec; 
            return *this; 
        }

        Ellipticity& operator=(Moments const & moments);
        Ellipticity& operator=(Axes const & axes);
    };

    class Axes {
        Eigen::Vector3d _data;
        // swap a,b and rotate if a<b, ensure theta in [-pi/2,pi/2)
        void normalize();     
    public:

        void setMajorAxis(double a) { 
            _data.x() = a; 
            normalize(); 
        }
        double getMajorAxis() const { return _data.x(); }

        void setMinorAxis(double b) { 
            _data.y() = b;
            normalize(); 
        }
        double getMinorAxis() const { return _data.y(); }

        void setAngle(double theta) { 
            _data.z() = theta; 
            normalize(); 
        }
        double getAngle() const { return _data.z(); }

        void set(double a, double b) { 
            _data.x() = a; 
            _data.y() = b; 
            normalize(); 
        }

        double getAxisRatio() const { return _data.y() / _data.x(); }

        template <typename Derived>
        explicit Axes(Eigen::MatrixBase<Derived> const & data) 
                : _data(data) { 
        }

        explicit Axes(double a=0, double b=0, double theta=0) 
                : _data(a,b,theta) {
        }

        Axes(Moments const & moments);
        Axes(Ellipticity const & ellipticity);

        template <typename Derived>
        void fill(Eigen::MatrixBase<Derived> & vector) const { vector = _data; }

        template <typename Derived>
        Axes& operator=(Eigen::MatrixBase<Derived> const & vec) { 
            _data = vec; 
            return *this; 
        }

        Axes& operator=(Moments const & moments);        
        Axes& operator=(Ellipticity const & ellipticity);

    };

    explicit Ellipse(Moments const & moments=Moments(), 
            Coordinate const & center=Coordinate(0,0))
            : _center(center), _moments(moments) {
    }

    explicit Ellipse(Axes const & axes, 
            Coordinate const & center=Coordinate(0,0))
            : _center(center), _moments(axes) {
    }
    
    explicit Ellipse(Ellipticity const & ellipticity, 
            Coordinate const & center=Coordinate(0,0)) 
        : _center(center), _moments(ellipticity) {        
    }

    Coordinate const & getCenter() const { return _center; }
    Coordinate & getCenter() { return _center; }
    
    Moments const & getMoments() const { return _moments; }
    Moments & getMoments() { return _moments; }
    void set(Moments const & moments) { _moments = moments; }

    Axes const getAxes() const { return Axes(_moments); }
    void set(Axes const & axes) { _moments = axes; }

    Ellipticity const getEllipticity() const { return Ellipticity(_moments); }
    void set(Ellipticity const & ellipticity) { _moments = ellipticity; }

    Ellipse transform(Eigen::Transform2d const & transform) const {
        return Ellipse(_moments, _center);
    }

    /** 
     *  \brief Functor that returns points on the boundary of the ellipse as a 
     *  function of a parameter that runs between 0 and 2pi (but is not angle).
     */
    class Parametric {
        Coordinate _center;
        Coordinate _u;
        Coordinate _v;
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Parametric(Ellipse const& ellipse);

        Coordinate operator()(double t) const {
            return _center + _u*std::cos(t) + _v*std::sin(t); 
        }
    };
    
    
    template <typename TRowCore, typename TColCore>
    struct ParameterJacobian {
        template <typename TMatrix>
        static void evaluate(
                const TColCore& core, 
                Eigen::MatrixBase<TMatrix>& m
        );
    };

    template <typename TRowCore>
    struct ParameterJacobian<TRowCore,TRowCore> {
        template <typename TMatrix>
        static void evaluate(
                const TRowCore& core, 
                Eigen::MatrixBase<TMatrix>& m
        ) {
        	m.setIdentity();
        }
    };



    template <typename TRowCore>
    struct ParameterJacobian<TRowCore,Moments> {
        template <typename TMatrix>
        static void evaluate(
                const Moments& core, 
                Eigen::MatrixBase<TMatrix>& m
        ) {
        	Eigen::Matrix3d tmp;
        	ParameterJacobian<Moments,TRowCore>::evaluate(core,tmp);
        	m = tmp.inverse();
        }
    };



private:
    Coordinate _center;
    Moments _moments;
};

template <>
struct Ellipse::ParameterJacobian<Ellipse::Moments, Ellipse::Moments> {
    template <typename TMatrix>
    static void evaluate(
            const Moments& core, 
            Eigen::MatrixBase<TMatrix>& m
    ) {
        m.setIdentity();
    }
};

template <>
struct Ellipse::ParameterJacobian<Ellipse::Moments, Ellipse::Axes> {
    template <typename TMatrix>
    static void evaluate(
            const Axes& core, 
            Eigen::MatrixBase<TMatrix>& m
    ) {
        double theta = core.getAngle();
        double a = core.getMajorAxis();
        double b = core.getMinorAxis();
        double c = std::cos(2.0*theta);
        double s = std::sin(2.0*theta);
        double cp = 1.0+c;
        double cm = 1.0-c;
        double a2 = a*a;
        double b2 = b*b;
        m(0,0) = a*cp;  m(0,1) = b*cm;  m(0,2) = (b2-a2)*s;
        m(1,0) = a*cm;  m(1,1) = b*cp;  m(1,2) = (a2-b2)*s;
        m(2,0) = a*s;   m(2,1) = -b*s;  m(2,2) = (a2-b2)*c;
    }
};

template <>
struct Ellipse::ParameterJacobian<Ellipse::Moments, Ellipse::Ellipticity> {
    template <typename TMatrix>
    static void evaluate(
            const Ellipticity& core, 
            Eigen::MatrixBase<TMatrix>& m
    ) {
        double e1 = core.getE1();
	    double e2 = core.getE2();
        double r = core.getRadius();
	    double norm = e1*e1 + e2*e2;
	    double e1p1 = e1 + 1;
	    double e1m1 = e1 - 1;
	    double r2 = r*r;
	    double n_inv = 1.0 / std::sqrt(1 - norm);
	    double n_inv_3 = n_inv*n_inv*n_inv;
        m(0,0) = (e1p1 - e2*e2)*r2*n_inv_3;  
        m(0,1) = e1p1*e2*r2*n_inv_3;     
        m(0,2) = 2*e1p1*r*n_inv;
        m(1,0) = (e1m1 + e2*e2)*r2*n_inv_3;  
        m(1,1) = -e1m1*e2*r2*n_inv_3;
        m(1,2) = -2*e1m1*r*n_inv;
        m(2,0) = e1*e2*r2*n_inv_3;
        m(2,1) = -e1m1*e1p1*r2*n_inv_3;
        m(2,2) = 2*e2*r*n_inv;
    }
};

template <typename TRowCore, typename TColCore>
template <typename TMatrix>
inline void Ellipse::ParameterJacobian<TRowCore, TColCore>::evaluate(
        const TColCore& core, 
        Eigen::MatrixBase<TMatrix>& m
) {
    Eigen::Matrix3d tmp;
    ParameterJacobian<TRowCore,Moments>::evaluate(core,m);
    ParameterJacobian<Moments,TColCore>::evaluate(core,tmp);
    m *= tmp;
}

}}} //end namespace lsst::meas::multifit

#endif // LSST_MEAS_MULTIFIT_ELLIPSE_H
