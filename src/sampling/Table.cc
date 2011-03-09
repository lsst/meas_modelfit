#include "lsst/meas/multifit/sampling/Table.h"

namespace ndd = lsst::ndarray::detail;

namespace lsst { namespace meas { namespace multifit { namespace sampling {

Table Table::allocate(
    int tableSize, int parameterSize, NestedSampleType nestedSampleType, int nestedSize
) {
    int stride = 3 + parameterSize;
    if (nestedSize > 0) stride += 1 + (nestedSize + 1) * nestedSize;
    std::pair<ndarray::Manager::Ptr,double*> r = ndarray::SimpleManager<double>::allocate(tableSize * stride);
    ndd::Core<1>::ConstPtr coreT = ndd::Core<1>::create(
        ndarray::makeVector(tableSize),
        ndarray::makeVector(stride),
        r.first
    );
    ndd::Core<2>::ConstPtr coreTP = ndd::Core<2>::create(
        ndarray::makeVector(tableSize, parameterSize), 
        ndarray::makeVector(stride, 1),
        r.first
    );
    ndd::Core<2>::ConstPtr coreTN = ndd::Core<2>::create(
        ndarray::makeVector(tableSize, nestedSize), 
        ndarray::makeVector(stride, 1),
        r.first
    );
    ndd::Core<3>::ConstPtr coreTNN = ndd::Core<3>::create(
        ndarray::makeVector(tableSize, nestedSize, nestedSize), 
        ndarray::makeVector(stride, nestedSize, 1),
        r.first
    );
    // if we don't have a nested expansion, set &nested.scalar == &target.
    double * p = r.second + ((nestedSize > 0) ? (3 + parameterSize) : 1);
    return Table(
        ndd::ArrayAccess< ndarray::Array<double,1,0> >::construct(r.second, coreT),
        ndd::ArrayAccess< ndarray::Array<double,1,0> >::construct(r.second + 1, coreT),
        ndd::ArrayAccess< ndarray::Array<double,1,0> >::construct(r.second + 2, coreT),
        ndd::ArrayAccess< ndarray::Array<double,2,1> >::construct(r.second + 3, coreTP),
        ndd::ArrayAccess< ndarray::Array<double,1,0> >::construct(p, coreT),
        ndd::ArrayAccess< ndarray::Array<double,2,1> >::construct(p + 1, coreTN),
        ndd::ArrayAccess< ndarray::Array<double,3,2> >::construct(p + nestedSize + 1, coreTNN),
        0.0, nestedSampleType
    );
}

void Table::clear() const {
    this->proposal = 0.0;
    this->target = 0.0;
    this->weight = 0.0;
    this->parameters = 0.0;
    this->_weightSum = 0.0;
    if (getNestedSize() > 0) {
        this->nested.scalar = 0.0;
        this->nested.vector = 0.0;
        this->nested.matrix = 0.0;
    }
}

void Table::computeWeights() const {
    weight = target;
    weight /= proposal;
    _weightSum = ndarray::sum(weight);
}

template class detail::RecordBase<double>;
template class detail::RecordBase<double const>;

template class detail::TableBase<double>;
template class detail::TableBase<double const>;

}}}} // namespace lsst::meas::multifit::sampling
