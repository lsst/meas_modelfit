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

#ifndef LSST_MEAS_MULTIFIT_Evaluator
#define LSST_MEAS_MULTIFIT_Evaluator

#include "lsst/base.h"
#include "lsst/meas/multifit/Grid.h"
#include "lsst/meas/multifit/BaseEvaluator.h"

namespace lsst { namespace meas { namespace multifit {

class Grid;

class Evaluator : public BaseEvaluator {
public:

    typedef boost::shared_ptr<Evaluator> Ptr;

    Wcs::Ptr getWCS() const {return _grid->wcs;}
#ifndef SWIG
    Grid::ConstPtr getGrid() const {return _grid;}

    Definition makeDefinition() const;
    Definition makeDefinition(
        ndarray::Array<double const,1,1> const & parameters
    ) const;
    
    static Ptr make(Definition const & definition);
#endif
   
    template<typename PixelT>
    static Ptr make(
        PTR(afw::image::Exposure<PixelT>) const & exposure,
        Footprint::Ptr const & fp,
        afw::geom::Point2D const & position,
        bool isVariable=false,
        bool fixPosition=true
    );

    template<typename PixelT>
    static Ptr make(
        PTR(afw::image::Exposure<PixelT>) const & exposure,
        Footprint::Ptr const & fp,
        ModelBasis::Ptr const & basis,
        afw::geom::ellipses::Ellipse const & ellipse,
        bool fixEllipticity=true,
        bool fixRadius=true,
        bool fixPosition=true
    );


protected:

    virtual void _evaluateModelMatrix(
        ndarray::Array<double,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & param
    ) const;

    virtual void _evaluateModelDerivative(
        ndarray::Array<double,3,3> const & matrix,
        ndarray::Array<double const,1,1> const & param
    ) const;

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & param) const;

private:
    
    FRIEND_MAKE_SHARED_1(Evaluator, boost::shared_ptr<lsst::meas::multifit::Grid>);

    explicit Evaluator(Grid::Ptr const & grid);

    Evaluator(Evaluator const & other);
    
    void _initialize();

    Grid::Ptr _grid;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluator
