#ifndef LSST_MEAS_MULTIFIT_SERSIC_CACHE_H
#define LSST_MEAS_MULTIFIT_SERSIC_CACHE_H


#include <Eigen/Core>
#include "lsst/afw/geom.h"
#include "lsst/pex/policy/DefaultPolicyFile.h"

namespace boost {
namespace serialization {
    class access;
}}

namespace lsst {
namespace meas {
namespace multifit {

class InterpolationFunction {
public:
    enum ExtrapolationParams {A=0, B, C, NPARAM};

    typedef boost::shared_ptr<InterpolationFunction> Ptr;
    typedef boost::shared_ptr<InterpolationFunction const> ConstPtr;
    typedef Eigen::Matrix<double, 1, NPARAM> ExtrapolationParameters;

    InterpolationFunction(
        double const min, double const max, double const step,
        Eigen::VectorXd const & values,
        ExtrapolationParameters const & extrapolationParameters
    ) : _min(min), _max(max), _step(step), 
        _values(values), 
        _extraParams(extrapolationParameters) 
    {}

    double operator() (double x) const;
    double d(double x) const;

protected:
    double _min, _max, _step;
    Eigen::VectorXd _values;
    ExtrapolationParameters _extraParams;
};

class SersicCache {
public:   
    typedef boost::shared_ptr<SersicCache> Ptr;
    typedef boost::shared_ptr<SersicCache const> ConstPtr;

    typedef InterpolationFunction Interpolator;

    static ConstPtr get(std::string const & name="");
    static ConstPtr make(
        lsst::pex::policy::Policy policy, 
        bool const doOverwrite=false,
        std::string const & name=""
    );
    static ConstPtr make(
        double const sersicMin, double const sersicMax,
        double const kMin, double const kMax,
        int const nSersic, int const nK,        
        double const innerRadius, double const outerRadius,
        double const extrapolationK1, double const extrapolationK2,
        double epsrel=-1.0, double epsabs=-1.0,        
        bool const & doOverwrite=false,
        std::string const & name="" 
    );
    static ConstPtr load(
        std::string const & filepath,         
        bool const & doOverwrite=false,
        std::string const & name="" 
    ); 
    void save(std::string const & filepath) const;

    Interpolator::ConstPtr getInterpolator(double const sersic) const;
    Interpolator::ConstPtr getDerivativeInterpolator(double const sersic) const;

    double const getSersicMin() const {return _sersicMin;}
    double const getSersicMax() const {return _sersicMax;}
    double const getKMin() const {return _kMin;}
    double const getKMax() const {return _kMax;}
    double const getSersicStep() const {return _sersicStep;}
    double const getKStep() const {return _kStep;}
    double const getInnerSersicRadius() const {return _inner;}
    double const getOuterSersicRadius() const {return _outer;}
    double const getExtrapolationK1() const {return _extrapolationK1;}
    double const getExtrapolationK2() const {return _extrapolationK2;}
    double const getEpsrel() const {return _epsrel;}
    double const getEpsabs() const {return _epsabs;}
    int const getNSersic() const {return int((_sersicMax - _sersicMin)/_sersicStep);}
    int const getNK() const {return int((_kMax - _kMin)/_kStep);}
    
    Eigen::MatrixXd const & getDataPoints() const {return _dataPoints;}
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & ar, unsigned int const);
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Interpolator::NPARAM> ExtrapolationParameters;

    SersicCache(
        double const sersicMin, double const sersicMax,
        double const kMin, double const kMax,
        int const nSersic, int const nK,        
        double const innerRadius, double const outerRadius,
        double const extrapolationK1, double const extrapolationK2,
        double epsrel=-1.0, double epsabs=-1.0
    );

    SersicCache(){};



    double _sersicMin, _sersicMax, _kMin, _kMax;
    double _sersicStep, _kStep;    
    double _inner, _outer;
    double _extrapolationK1, _extrapolationK2;
    double _epsrel, _epsabs;
    Eigen::MatrixXd _dataPoints;
    ExtrapolationParameters _extraParams;
    static std::map<std::string, Ptr> _registry;
};
}}}
#endif


