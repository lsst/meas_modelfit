#ifndef LSST_MEAS_MULTIFIT_MODELITERATOR_H
#define LSST_MEAS_MULTIFIT_MODELITERATOR_H

#include "boost/iterator/iterator_facade.hpp"

#include "lsst/meas/multifit/Model.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 *  Specialized iterator over Model::Map.  Dereferences to Model const &,
 *  and supplies getTag(), getLinearOffset(), and getNonlinearOffset()
 *  to retreive the associated tag and parameter indices into a multi-model
 *  vector or matrix.
 */
class ObjectModelIterator : public boost::iterator_facade<
        ObjectModelIterator,
        ObjectModel::Ptr,
        boost::bidirectional_traversal_tag,
        ObjectModel::Ptr>
{
    typedef ObjectModel::Map::iterator BaseIterator;

    BaseIterator _base;
    int _linear_offset;
    int _nonlinear_offset;

    friend class boost::iterator_core_access;

    ObjectModel::Ptr dereference() const { return _base->second; }

    void increment() {
        ObjectModel::Ptr model = dereference();
        _linear_offset += model->getNumLinearParam();
        _nonlinear_offset += model->getNumNonlinearParam();
        ++_base;
    }

    void decrement() {
        --_base;
        ObjectModel::Ptr model = dereference();
        _linear_offset -= model->getNumLinearParam();
        _nonlinear_offset -= model->getNumNonlinearParam();
    }
    bool equal(ObjectModelIterator const & other) const { 
        return _base == other._base; 
    }

public:

    ObjectModelIterator()
        : _base(), _linear_offset(0), _nonlinear_offset(0) {}

    explicit ObjectModelIterator(BaseIterator const & base) 
        : _base(base), _linear_offset(0), _nonlinear_offset(0) {}

    ObjectModelIterator(ObjectModelIterator const & other) 
        : _base(other._base), _linear_offset(other._linear_offset), 
          _nonlinear_offset(other._nonlinear_offset) {}

    int getLinearOffset() const { return _linear_offset; }
    int getNonlinearOffset() const { return _nonlinear_offset; }

    const ObjectModel::Key & getObjectId() const { return _base->first; }

    bool operator==(BaseIterator const & other) const { 
        return this->_base == other; 
    }
    bool operator!=(BaseIterator const & other) const { 
        return this->_base != other; 
    }
    static ObjectModelIterator begin(ObjectModel::Map & models) { 
        return ObjectModelIterator(models.begin()); 
    }
    static ObjectModelIterator end(ObjectModel::Map & models) {
        return ObjectModelIterator(models.end());    
    }
    static ObjectModelIterator find(ObjectModel::Map & models, 
            ObjectModel::Key const objectId) {
        ObjectModelIterator r = begin(models);
        while (r != models.end()) {
            if (r.getObjectId() == objectId) return r;
            ++r;
        }
        return r;
    }
};

}}} //end namespace lsst::meas::multifit

#endif
