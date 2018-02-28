#include <vector>
#include <Eigen/Dense>

namespace lsst {
namespace meas {
namespace modelfit {

struct Moments {
    using FirstMoment = Eigen::Matrix<double, 2, 1>
    using SecondMoment = Eigen::Matrix<double, 2, 2>
    using ParameterVector = Eigen::Matrix<double, 6, 1>
    using ParameterMatrix = Eigen::Matrix<double, 6, 6>

    template <typename iter>
    Moments(iter begin, iter end);

    Moments(std::vector<double> moments);

    Moments(double zero, FirstMoment first, SecondMoment second);

    double operator[](int i);

    ParameterVector getParameterVector();

private:
    double zero;
    FirstMoment one;
    SecondMoment two;
}

struct ZerothShapeMoment {
    static double computeValue(Moments Q, Moments W);

    static Moments::ParameterVector computeGradient(Moments Q, Moments W);
};

struct FirstShapeMoment {
    static double computeValue(Moments Q, Moments W);

    static Moments::ParameterVector computeGradient(Moments Q, Moments W);
};

struct SecondShapeMoment {
    static double computeValue(Moments Q, Moments W);

    static Moments::ParameterVector computeGradient(Moments Q, Moments W);
};
}}} // Close namespace lsst::meas::modelfit
