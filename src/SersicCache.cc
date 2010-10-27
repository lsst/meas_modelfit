#include "lsst/meas/multifit/SersicCache.h"
#include "lsst/meas/multifit/ModifiedSersic.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Debug.h"
#include "lsst/utils/ieee.h"

#include <fstream>
#include <boost/serialization/list.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/make_shared.hpp>
#include <Eigen/LU>
#include <Eigen/Array>

namespace multifit = lsst::meas::multifit;

namespace {

void fitRational(
    Eigen::Vector3d const & k, 
    Eigen::Vector3d const & values,
    double & a0, double & a1, double & a2
) {
    Eigen::Matrix3d matrix;
    matrix.col(0).fill(1.0);
    matrix.col(1) = k.cwise().inverse();
    matrix.col(2) = k.cwise().inverse().cwise().square();
    Eigen::Vector3d model;
    Eigen::LU<Eigen::Matrix3d> lu(matrix);
    lu.solve(values, &model);
    a0 = model[0];
    a1 = model[1];
    a2 = model[2];
    if (lsst::utils::isnan(a0) || lsst::utils::isnan(a1) || lsst::utils::isnan(a2)) {
        a0 = a1 = a2 = 0.0;        
    }
}

} //end annonymous namespace


///////////////////////////////////////////////////////////////////////////////
// InterpolationFunction
///////////////////////////////////////////////////////////////////////////////
double multifit::InterpolationFunction::operator()(double p) const {
    if (p < _max) {
        if (p < _min) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                    p % _min % _max).str()
            );
        }
        int i = static_cast<int>((p - _min) / _step);    


        double w1 = (p - (_min + _step*i))/_step;
        double w0 = 1 - w1;
        double v0 = _values[i];
        double v1 = _values[i+1];

        return w0*v0 + w1*v1;
    }
    
    return _extraParams[A] + _extraParams[B] / p + _extraParams[C] / (p*p);
}

double multifit::InterpolationFunction::d(double p) const {
    if (p < _max) {
       if (p < _min) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                (boost::format("Operand value (%1%) is outside of valid range [%2%, %3%)")%
                    p % _min % _max).str()
            );
        }
        int i = static_cast<int>((p - _min) / _step);    

        double v0 = _values[i];
        double v1 = _values[i+1];
         
        return (v1 - v0)/_step; 
    }
    return -_extraParams[B] / (p*p) - 2.0 * _extraParams[C] / (p*p*p);
}
///////////////////////////////////////////////////////////////////////////////
// SersicCache
///////////////////////////////////////////////////////////////////////////////

std::map<std::string, multifit::SersicCache::Ptr> multifit::SersicCache::_registry;

multifit::SersicCache::ConstPtr multifit::SersicCache::get(std::string const & name) {
    return _registry[name];
}

multifit::SersicCache::ConstPtr multifit::SersicCache::make(
    lsst::pex::policy::Policy policy,
    bool const doOverwrite,
    std::string const & name
) {
    lsst::pex::policy::DefaultPolicyFile defSource(
        "meas_multifit",
        "SersicCacheDict.paf",
        "policy"
    );
    lsst::pex::policy::Policy defPol(defSource);
    policy.mergeDefaults(defPol);

    return SersicCache::make(
        policy.getDouble("sersicMin"), policy.getDouble("sersicMax"),
        policy.getDouble("kMin"), policy.getDouble("kMax"),
        policy.getInt("nBinSersic"), policy.getInt("nBinK"),
        policy.getDouble("sersicInnerRadius"), 
        policy.getDouble("sersicOuterRadius"), 
        policy.getDouble("extrapolationK1"),
        policy.getDouble("extrapolationK2"),
        policy.getDouble("epsabs"), 
        policy.getDouble("epsrel"),
        doOverwrite,
        name
    );
}


multifit::SersicCache::ConstPtr multifit::SersicCache::make(
    double const sersicMin, double const sersicMax,
    double const kMin, double const kMax,
    int const nSersic, int const nK,
    double const innerRadius, double const outerRadius,
    double const extrapolationk1, double const extrapolationk2,
    double epsrel, double epsabs,
    bool const & doOverwrite,
    std::string const & name
) {

    if(!doOverwrite && !name.empty()) {
        std::map<std::string, Ptr>::const_iterator i = _registry.find(name);
        if(i != _registry.end()){
            return i->second;
        }
    }
    Ptr cache(
        new SersicCache(
            sersicMin, sersicMax, kMin, kMax, 
            nSersic, nK, 
            innerRadius, outerRadius,
            extrapolationk1, extrapolationk2,
            epsrel, epsabs
        )
    );
    if(doOverwrite && !name.empty()) {
            _registry[name] = cache;
    }
    
    return cache;
}

multifit::SersicCache::ConstPtr multifit::SersicCache::load(
    std::string const & filepath,
    bool const & doOverwrite,
    std::string const & name
) {
    std::map<std::string, Ptr>::const_iterator i = _registry.find(name);
    bool useRegistry = i != _registry.end() && !doOverwrite;
    if(useRegistry) {
        return i->second;        
    }    

    Ptr cache(new SersicCache());
    std::ifstream ifs(filepath.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> *cache;
    
    if(!useRegistry) {
        _registry[name] = cache;
    }

    return cache;
}


void multifit::SersicCache::save(std::string const & filepath) const {
    std::ofstream ofs(filepath.c_str());
    boost::archive::text_oarchive oa(ofs);
    oa << *this;
}

multifit::SersicCache::SersicCache(
    double const sersicMin, double const sersicMax,
    double const kMin, double const kMax,
    int const nSersic, int const nK,
    double const innerRadius, double const outerRadius,
    double const extrapolationK1, double const extrapolationK2,
    double epsrel, double epsabs
) : _sersicMin(sersicMin), _sersicMax(sersicMax), _kMin(kMin), _kMax(kMax),
    _sersicStep((sersicMax-sersicMin)/nSersic), _kStep((kMax-kMin)/nK),
    _inner(innerRadius), _outer(outerRadius),
    _extrapolationK1(extrapolationK1), _extrapolationK2(extrapolationK2),
    _epsrel(epsrel), _epsabs(epsabs),
    _dataPoints(nSersic+1, nK+1),
    _extraParams(nSersic+1, Interpolator::NPARAM)
{
    lsst::pex::logging::Debug debug("lsst.meas.multifit.SersicCache");

    //compute the data points
    //outer loop over sersic dimension
    int nRows = nSersic+1, nCols = nK+1;

    double sersic=sersicMin;
    ModifiedSersicHankelTransform hankel(
        ModifiedSersicFunction(sersicMin, _inner, _outer),
        epsrel, epsabs
    );
    for(int i = 0; i < nRows; ++i, sersic+=_sersicStep) {
        if(i > 0) {
            hankel.setSersicIndex(sersic);
        }
        double k = kMin;
        debug.debug<5>("Filling row %d of %d.", i, nRows);
        //inner loop over k dimension 
        Eigen::MatrixXd::RowXpr row = _dataPoints.row(i);
        for (int j = 0; j < nCols; ++j, k += _kStep) {
            //set grid point
            debug.debug<8>("Filling item %d of %d (column %d of %d).",
                           i*nCols + j, nCols * nRows, j, nCols);
            row[j] = hankel(k);
        }
        
        double v0 = row[nK];
        double v1 = hankel(extrapolationK1);
        double v2 = hankel(extrapolationK2);
        
        ExtrapolationParameters::RowXpr params = _extraParams.row(i);
        ::fitRational(
            Eigen::Vector3d(kMax, extrapolationK1, extrapolationK2),
            Eigen::Vector3d(v0, v1, v2),
            params[Interpolator::A], params[Interpolator::B], params[Interpolator::C]
        );
    }
}

multifit::SersicCache::Interpolator::ConstPtr multifit::SersicCache::getInterpolator(
    double const sersic
) const {
    if (sersic < _sersicMin || sersic >= _sersicMax) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                sersic % _sersicMin % _sersicMax).str()
        );
    }

    int i = static_cast<int>((sersic - _sersicMin) / _sersicStep);
    double s = (sersic - (_sersicMin+_sersicStep*i))/_sersicStep;
        
    Eigen::VectorXd r(_dataPoints.row(i)*(1 - s) + _dataPoints.row(i+1)*s);
    Interpolator::ExtrapolationParameters params =
        _extraParams.row(i)*(1-s) + _extraParams.row(i+1)*s;
    return boost::make_shared<Interpolator>(_kMin, _kMax, _kStep, r, params);
}

multifit::SersicCache::Interpolator::ConstPtr multifit::SersicCache::getDerivativeInterpolator(
    double const sersic
) const {
    if (sersic < _sersicMin || sersic >= _sersicMax) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Operand value %1% is outside of valid range [%2%, %3%)")%
                sersic % _sersicMin % _sersicMax).str()
        );
    }

    int i = static_cast<int>((sersic - _sersicMin) / _sersicStep);
    Eigen::VectorXd r((_dataPoints.row(i+1) - _dataPoints.row(i))/ _sersicStep);
    Interpolator::ExtrapolationParameters params =
        (_extraParams.row(i)*(i+1) + _extraParams.row(i))/_sersicStep;

    return boost::make_shared<Interpolator>(_kMin, _kMax, _kStep, r, params);
}

template <class Archive>
void multifit::SersicCache::serialize(Archive & ar, unsigned int const) {
    ar & _sersicStep;
    ar & _kStep;
    ar & _sersicMin;
    ar & _sersicMax;
    ar & _kMin;
    ar & _kMax;
    ar & _inner;
    ar & _outer;
    ar & _epsrel;
    ar & _epsabs;
    ar & _extrapolationK1;
    ar & _extrapolationK2;

    if (Archive::is_loading::value) {
        int nRow = getNSersic() + 1;
        int nCol = getNK() + 1;

        _dataPoints = Eigen::MatrixXd(nRow, nCol);
        _extraParams = ExtrapolationParameters(nRow, Interpolator::NPARAM);
    }
    
    double * iData = _dataPoints.data();
    for(int i =0; i < _dataPoints.size(); ++i, ++iData) {
        ar & *iData;
    }
    double * iParam = _extraParams.data();
    for(int i =0; i < _extraParams.size(); ++i, ++iParam) {
        ar & *iParam;
    }
}

template void multifit::SersicCache::serialize<boost::archive::text_oarchive> (
    boost::archive::text_oarchive &, unsigned int const
);
template void multifit::SersicCache::serialize<boost::archive::text_iarchive> (
    boost::archive::text_iarchive &, unsigned int const
);


