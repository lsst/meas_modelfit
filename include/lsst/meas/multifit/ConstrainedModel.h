#ifndef LSST_MEAS_MULTIFIT_CONSTRAINEDMODEL_H
#define LSST_MEAS_MULTIFIT_CONSTRAINEDMODEL_H

#include "ndarray/ndarray.hpp"
#include "lsst/meas/multifit/Model.h"

namespace lsst {
namespace meas {
namespace multifit {

// constrained model that allows all linear parameters to vary together.
// NOTE: this template does not register the constrained model with the factory;
// this must be done explicitly for each model that is to be constrained.
template <typename TObjectModel>
class FluxConstrainedModel : public TObjectModel {

public:
    typedef TObjectModel Base;
    typedef ObjectModel::SourceModel SourceModel;
    
    virtual int getNumLinearParam() const { return 1; }

    MULTIFIT_MODEL_CONSTRAINT("flux")

    FluxConstrainedModel(
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinear_params, 
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & linear_params)
    : Base(nonlinear_params,linear_params), _tied_linear_params(linear_params) {    }

    virtual void computeLinearMatrix(SourceModel::Ptr source,
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinear_params,
            ndarray::ArrayRef<ObjectModel::Real,3,2> const & matrix_linear,
            ndarray::ArrayRef<ObjectModel::Real,2,2> const & k) const {
        assert(matrix_linear.shape()[0] == 1);
        ndarray::Array<ObjectModel::Real,3> full_matrix_linear(
                ndarray::make_index(_tied_linear_params.size(),
                        matrix_linear.shape()[1],
                        matrix_linear.shape()[2])
        );
        Base::computeLinearMatrix(source,nonlinear_params,full_matrix_linear,k);
        //todo fix implementation below
        #if 0
        matrix_linear.core().matrix() += full_matrix_linear.core().matrix(
                this->_tied_linear_params.size(),
                matrix_linear.shape()[1]*-_tied_linear_params.shape()[2]).transpose()
                * this->_tied_linear_params.core().vector();
        #endif
    }

    virtual void computeNonlinearMatrix(SourceModel& source,
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & linear_params,
            ndarray::ArrayRef<ObjectModel::Real,3,2> const & m_nonlinear) const 
    {
        assert(linear_params.size() == 1);
        ndarray::Array<ObjectModel::Real,1> full_linear_params(
                _tied_linear_params.size()
        );
        full_linear_params.core().vector() = 
                linear_params[0] * this->_tied_linear_params.core().vector();
        Base::computeNonlinearMatrix(source,full_linear_params,m_nonlinear);
    }

private:
    ndarray::Array<ObjectModel::Real,1> _tied_linear_params;
};

}}} //end namespace lsst::meas::multifit

#endif
