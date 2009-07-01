#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include <map>
#include <stdexcept>

#include "boost/shared_ptr.hpp"
#include "boost/utility.hpp"
#include "boost/cstdint.hpp"

#include "Eigen/Core"
#include "ndarray/ndarray.hpp"

#include "lsst/meas/multifit/CalibratedExposure.h"

namespace Eigen {
    typedef Matrix<double, 5, 1> Vector5d;    
}

namespace lsst {
namespace meas {
namespace multifit {

#define MULTIFIT_MODEL_TYPE(type)                                       \
    static std::string getTypeNameStatic() { return type; }             \
    virtual std::string getTypeName() const { return getTypeNameStatic(); }

#define MULTIFIT_MODEL_CONSTRAINT(constraint)                           \
    static std::string getConstraintNameStatic() { return constraint; } \
    virtual std::string getConstraintName() const {                     \
        return getConstraintNameStatic();                               \
    }

#define MULTIFIT_MODEL_REGISTER(name)                                   \
    static ModelFactory<name>                                           \
    name ## _factory(name::getTypeNameStatic(),name::getConstraintNameStatic())



/** 
 *  Subclasses of ObjectModel define how vectors of nonlinear and linear 
 *  parameters (which are not owned by the model) are interpreted onto a 
 *  pixelized model image, and possibly how they are to be constrained during 
 *  fitting. ObjectModel and its subclasses should be completely immutable, 
 *  and in general will not have much data (often only the expected size of its 
 *  parameter vectors, when this is not known statically).
 *
 *  Special considerations concrete subclasseses (here called ModelT):
 *
 *   - A constructor with the same signature as ModelFactoryBase::create must be
 *     provided.Here the parameter vectors are the full, unconstrained parameter
 *     vectors, and hence may be different from the constrained parameter 
 *     vectors that SourceModels will receive.
 *
 *   - A static member or static global variable of type ModelFactory<ModelT> 
 *     must be declared an initialized with the type and constraint name values 
 *     corresponding to the model. This can be handled by using the 
 *     MULTIFIT_MODEL_REGISTER macro.
 *
 *   - A concrete subclass of ObjectModel::SourceModel must be present as 
 *     ModelT::SourceModel. This may be a typedef or nested class, but it 
 *     should only be constructable via createSource and a protected 
 *     constructor that takes only a single Exposure parameter (and 
 *     createSource should call that constructor).
 *
 *   - The getConstraintName() and getTypeName() methods should refer to static
 *     member functions getConstraintNameStatic() and getTypeNameStatic(), 
 *     respectively.  These can all be defined using the MULTIFIT_MODEL_TYPE 
 *     and MULTIFIT_MODEL_CONSTRAINT macros.
 */
class ObjectModel : boost::noncopyable {
public:
    typedef boost::shared_ptr<ObjectModel> Ptr;

    typedef double Real;
    typedef boost::int64_t Key;

    typedef std::map<Key const, ObjectModel::Ptr> Map;

    ndarray::ArrayRef<Real,1,1> extractNonlinearParameters(
            const ndarray::ArrayRef<Real,1,1> & allParameters) const {
        return allParameters[ndarray::Range(0,getNumNonlinearParam())];
    }

    ndarray::ArrayRef<Real,1,1> extractLinearParameters(
            const ndarray::ArrayRef<Real,1,1> & allParameters) const {
        int nl = getNumNonlinearParam();
        return allParameters[ndarray::Range(nl,nl+getNumLinearParam())];
    }

    virtual int getNumLinearParam() const = 0; 
    // This value must not change after construction.

    virtual int getNumNonlinearParam() const = 0; 
    // This value must not change after construction.

    virtual std::string getConstraintName() const = 0;
    virtual std::string getTypeName() const = 0;

    // SourceModels manage a workspace with temporaries related to evaluating 
    // the modelon a particular exposure. It is expected that 
    // computeLinearMatrix, computeNonlinearMatrix, and computeCalibrationMatrix    
    // will be called in that order.

    // ObjectModel derived classes must also define their own SourceModel, 
    // derived from ObjectModel::SourceModel. This definition may be a typedef 
    // to a non-abstract base class SourceModel.
    class SourceModel : boost::noncopyable {

    public:
        typedef boost::shared_ptr<SourceModel> Ptr;
        typedef std::pair<CalibratedExposure::Ptr const,ObjectModel::Key const> Key;
        typedef std::map<Key, SourceModel::Ptr> Map;

        const CalibratedExposure& getExposure() const { return *_exposure; }

        virtual ~SourceModel() {}

    protected:
        // constructor for all derived classes must have the same signature
        CalibratedExposure::Ptr _exposure;
        
        explicit SourceModel(
                CalibratedExposure::Ptr const & exposure,
                int marginalizaionFlags) 
            : _exposure(exposure) {}

        friend class ObjectModel;
    };

    // An override of this method should be the only way to create a 
    // SourceModel-derived class. It should not use any ObjectModel data to 
    // create the SourceModel (in other words, overrides should look exactly 
    // the same, only with "SourceModel" referring to a class derived from 
    // ObjectModel::SourceModel).
    virtual SourceModel* createSource(
            CalibratedExposure::Ptr const & exposure,
            int marginalizationFlags) const {
        return new SourceModel(exposure, marginalizationFlags);
    }

    // Calling computeLinearMatrix sets the nonlinear parameters to be used
    // by the SourceModel for later calls to computeNonlinearMatrix and 
    // computeCalibrationMatrix.  If nonlinear_parameters changes, 
    // computeLinearMatrix must be called again.
    virtual void computeLinearMatrix(SourceModel& source,
            ndarray::ArrayRef<Real,1,1> const & nonlinearParameters,
            ndarray::ArrayRef<Real,3,2> const & linearMatrix,
            ndarray::ArrayRef<Real,2,2> const & residuals) const = 0;
    
    // Calling computeNonlinearMatrix sets the linear parameters to be used
    // by the SourceModel for later calls to computeCalibrationMatrix.  
    // If linear_parameters changes, computeNonlinearMatrix must be called 
    // again.
    virtual void computeNonlinearMatrix(SourceModel& source, 
            ndarray::ArrayRef<Real,1,1> const & linearParameters,
            ndarray::ArrayRef<Real,3,2> const & nonlinearMatrix) const = 0;
    
    virtual void computePsfMatrix(SourceModel * source, 
            ndarray::ArrayRef<Real,3,2> const & psfMatrix) const = 0;

    virtual ~ObjectModel() {}

};

// Factory function for ObjectModel. All ObjectModel-derived classes should 
// be constructable from a type string, a constraint string, and a pair of 
// unconstrained parameter vectors.
// Concrete ObjectModels are registered with the factory simply by defining a 
// static ModelFactory member templated on that ObjectModel class. 
ObjectModel* createObjectModel(std::string const & type, 
        std::string const & constraint,
        ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinearParams,
        ndarray::ArrayRef<ObjectModel::Real,1,1> const & linearParams);

struct ModelFactoryBase {
public:
    virtual ~ModelFactoryBase() {}
    virtual ObjectModel * create(
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinearParams,
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & linearParams) 
            const = 0;
protected:
    ModelFactoryBase(std::string const & type, std::string const & constraint);
    
    typedef std::pair<std::string,std::string> RegistryKey;
    typedef std::map<RegistryKey, ModelFactoryBase *> RegistryMap; 
    typedef RegistryMap::iterator RegistryIterator;
    
    // saves us from static initialization order problems.
    static RegistryMap& getRegistry(); 

    friend ObjectModel * createObjectModel(
            std::string const & type, std::string const & constraint,
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinearParams, 
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & linearParams);
};

template<typename ModelT>
struct ModelFactory : public ModelFactoryBase {

    ModelT * create(
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinearParams, 
            ndarray::ArrayRef<ObjectModel::Real,1,1> const & linearParams) 
            const {
        return new ModelT(nonlinearParams,linearParams);
    }
    
    ModelFactory(std::string const & type, std::string const & constraint)
        : ModelFactoryBase(type,constraint) {}
    
};

}}} //end namespace lsst::meas::multifit
#endif
