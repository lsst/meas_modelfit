#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FIXED_MORPHOLOGY
#include "lsst/meas/multifit/components/Morphology.h"

template <typename BaseMorphology>
class FixedNonlinearMorphology : public BaseMorphology {
    typedef BaseMorphology Base;
    typedef boost::shared_ptr<Base> BasePtr;
    typedef boost::shared_ptr<Base const> BaseConstPtr;
    typedef FixedNonlinearMorphology<Base> Fixed;

    typedef boost::shared_ptr<ParameterVector> ParameterVectorPtr;
    typedef boost::shared_ptr<ParameterVector const> ParameterVectorConstPtr;
public:
    typedef boost::shared_ptr<Fixed> Ptr;
    typedef boost::shared_ptr<Fixed const> ConstPtr;

    static Ptr create(Base const & base) {
        ParameterVectorPtr nonlinearParameters;
        if (base.getNonlinearParameterSize() != 0) {
            nonlinearParameters.reset(
                new ParameterVector(base._getNonlinearParameterSize())
            );
            std::copy(
                base.beginNonlinear(), 
                base.endNonlinear(), 
                nonlinearParameters->data()
            );            
        } else {
            nonlinearParameters.reset(new ParameterVector());
        }

        return boost::make_shared<Fixed>(base.getLinearParameters(), nonlinearParameters);
    }
protected:     
    virtual int const getNonlinearParameterSize() const {return 0;}

    virtual _handleNonlinearParameterChanged() {};

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const &,
        boost::shared_ptr<ParameterVector const> const &
        size_t const &,
    ) const {
        return Fixed::create(*this);
    }
    
    FixedNonlinearMorphology(
        ParameterVectorConstPtr const & linearParameters,
        ParameterVectorConstPtr const & nonlinearParameters
    ) : Base(linearParameters, nonlinearParameters) {}
};
