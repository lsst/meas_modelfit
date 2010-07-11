/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
#ifndef LSST_MEAS_MULTIFIT_SERSIC_CACHE_H
#define LSST_MEAS_MULTIFIT_SERSIC_CACHE_H

#include <boost/make_shared.hpp>
#include "lsst/afw/geom.h"
#include "lsst/meas/multifit/Cache.h"
#include "lsst/pex/policy/DefaultPolicyFile.h"

namespace lsst {
namespace meas {
namespace multifit {

class SersicCache : public Cache {
public:
    typedef boost::shared_ptr<SersicCache> Ptr;
    typedef boost::shared_ptr<SersicCache const> ConstPtr;

    static ConstPtr getInstance();

private:
    typedef lsst::pex::policy::Policy::Ptr PolicyPtr;
    typedef boost::shared_ptr<lsst::pex::policy::PolicySource> PolicySourcePtr;

    class FillFunction : public Cache::FillFunction {
    public:        
        FillFunction (double const & epsabs, double const & epsrel, int const & limit) 
          : Cache::FillFunction(0),
            _lastY(std::numeric_limits<double>::quiet_NaN()),
            _epsabs(epsabs), _epsrel(epsrel), _limit(limit)
        {}

        virtual double operator()(double x, double y) const;
        virtual Cache::FillFunction::Ptr clone() const {
            return boost::make_shared<FillFunction>(_epsabs, _epsrel, _limit);
        }

    private:
        class IntegralParameters {
        public:    
            double const & getK() const {return _k;}
            double const & getN() const {return _n;}
            double const & getKappa() const {return _kappa;}
            double const & getNorm() const {return _norm;}

            void setK(double const & k) {_k = std::abs(k);}
            void setN(double const & n);      
        private:
            double _k, _n, _kappa, _norm;
        };

        static double sersicFunction(double radius, void * parameters);
            
        double _lastY;
        double _epsabs, _epsrel;
        int _limit;
        mutable IntegralParameters _params;
    };

    SersicCache(
        lsst::afw::geom::BoxD const & parameterBounds, 
        lsst::afw::geom::Extent2D const& resolution, 
        Cache::FillFunction::Ptr const & fillFunction
    ) : Cache(parameterBounds, resolution, fillFunction) {}

    static PolicySourcePtr getDefaultPolicySource() {
        static const PolicySourcePtr source(
            new lsst::pex::policy::DefaultPolicyFile(
                "meas_multifit",
                "SersicCacheDict.paf",
                "policy"
            )
        );
        return source;
    }

    static Ptr _singleton;
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_SERSIC_CACHE_H




