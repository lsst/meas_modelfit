// -*- lsst-c++ -*-

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
 
%define samplingLib_DOCSTRING
"
Basic routines to talk to lsst::meas::multifit::sampling classes
"
%enddef


%feature("autodoc", "1");
%module(package="lsst.meas.multifit.sampling", docstring=samplingLib_DOCSTRING) samplingLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored
#pragma SWIG nowarn=401                 // nothin known about base class X
%{
#include <boost/random/mersenne_twister.hpp>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/LocalPsf.h"
#include "lsst/meas/multifit/sampling/RandomEngine.h"
#include "lsst/meas/multifit/sampling/Table.h"
#include "lsst/meas/multifit/sampling/MixtureDistribution.h"
#include "lsst/meas/multifit/sampling/IterativeImportanceSampler.h"
#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_MULTIFIT_SAMPLING_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "lsst/ndarray/python.h"
#include "lsst/ndarray/python/eigen.h"
%}

%inline %{
namespace boost {}
namespace lsst { 
    namespace afw {
        namespace image {}
        namespace detection {}
        namespace math {}
        namespace geom {}
    }
namespace meas { namespace multifit { namespace sampling {
}}} // namespace meas::multifit::sampling
} // namespace lsst

using namespace lsst;
 using namespace lsst::meas::multifit::sampling;
%}

/******************************************************************************/
%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%include "std_complex.i"
%include "std_vector.i"

%lsst_exceptions();

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL$"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    version_svn = lsst.utils.guessSvnVersion(HeadURL)

    try:
        import eups
    except ImportError:
        return version_svn
    else:
        try:
            version_eups = eups.getSetupVersion("meas_multifit")
        except AttributeError:
            return version_svn

    if version_eups == version_svn:
        return version_svn
    else:
        return "%s (setup: %s)" % (version_svn, version_eups)
%}

%include "lsst/ndarray/ndarray.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/math/shapelets/shapeletsLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/meas/multifit/multifitLib.i"

%declareNumPyConverters(Eigen::VectorXd);
%declareNumPyConverters(Eigen::MatrixXd);

%declareNumPyConverters(lsst::ndarray::Array<double,1,1>);
%declareNumPyConverters(lsst::ndarray::Array<double,2,2>);

%declareNumPyConverters(lsst::ndarray::Array<double const,1,1>);
%declareNumPyConverters(lsst::ndarray::Array<double const,2,2>);

/*****************************************************************************/

%inline %{

    PyArray_Descr * makeTableNumpyType(int parameterSize, int nestedSize) {
        PyObject * nestedList = PyList_New(3);
        PyList_SET_ITEM(nestedList, 0, Py_BuildValue("sO", "scalar", &PyFloat_Type));
        PyList_SET_ITEM(nestedList, 1, Py_BuildValue("sOi", "vector", &PyFloat_Type, nestedSize));
        PyList_SET_ITEM(nestedList, 2, 
                        Py_BuildValue("sO(ii)", "matrix", &PyFloat_Type, nestedSize, nestedSize));
        PyArray_Descr * nestedType = 0;
        if (PyArray_DescrConverter(nestedList, &nestedType) < 0) {
            std::cerr << "checkpoint 1\n";
            return 0;
        }
        if (nestedType == 0) {
            std::cerr << "checkpoint 2\n";
            return 0;
        }
        Py_DECREF(nestedList);
        PyObject * mainList = PyList_New(5);
        PyList_SET_ITEM(mainList, 0, Py_BuildValue("sO", "importance", &PyFloat_Type));
        PyList_SET_ITEM(mainList, 1, Py_BuildValue("sO", "target", &PyFloat_Type));
        PyList_SET_ITEM(mainList, 2, Py_BuildValue("sO", "weight", &PyFloat_Type));
        PyList_SET_ITEM(mainList, 3, Py_BuildValue("sOi", "parameters", &PyFloat_Type, parameterSize));
        PyList_SET_ITEM(mainList, 4, Py_BuildValue("sN", "nested", nestedType));
        PyArray_Descr * mainType = 0;
        if (PyArray_DescrConverter(mainList, &mainType) < 0) {
            std::cerr << "checkpoint 3\n";
            return 0;
        }
        Py_DECREF(mainList);
        return mainType;
    }

    PyObject * returnTable(Table const & other) {
        PyArray_Descr * dtype = makeTableNumpyType(other.getParameterSize(), other.getNestedSize());
        if (dtype == 0) {
            std::cerr << "checkpoint 4\n";
            return 0;
        }
        PyObject * owner = PyCObject_FromVoidPtr(
            new lsst::ndarray::Manager::Ptr(other.importance.getManager()), 
            lsst::ndarray::detail::destroyCObject
        );
        npy_intp dim = other.getTableSize();
        npy_intp stride = other.importance.getStride<0>() * sizeof(double);
        PyObject * result = PyArray_NewFromDescr(
            &PyArray_Type, dtype, 1, &dim, &stride, other.importance.getData(),
            NPY_ALIGNED | NPY_WRITEABLE, 0
        );
        if (result == 0) {
            std::cerr << "checkpoint 5\n";
            return 0;
        }
        PyArray_BASE(result) = owner;
        return result;
    }

    PyObject * allocateTable(
        int tableSize, int parameterSize,
        lsst::meas::multifit::sampling::NestedSampleType sampleType, 
        int nestedSize
    ) {
        lsst::meas::multifit::sampling::Table r(
            Table::allocate(tableSize, parameterSize, sampleType, nestedSize)
        );
        return returnTable(r);
    }

    struct SampleIterator {
        typedef std::list<Table>::const_iterator CIter;

        CIter current;
        CIter end;

        SampleIterator() {}

        SampleIterator(CIter const & current_, CIter const & end_) : current(current_), end(end_) {}

        PyObject * next() {
            if (current == end) {
                PyErr_SetNone(PyExc_StopIteration);
                return 0;
            }
            PyObject * r = returnTable(*current);
            ++current;
            return r;
        }
    };

%}

namespace lsst { namespace meas { namespace multifit { namespace sampling {

enum NestedSampleType {
    MACLAURIN_SERIES, // ln(f(0)), [d ln(f) / dx](0), [d^2 ln(f) / dx^2](0)
    FISHER_FULL,      // ln(f(x_0)), x_0, [d^2 ln(f) / dx^2](x_0)
    FISHER_LLT,       // like FISHER_FULL, but stores the lower Cholesky factor of the matrix
    COVARIANCE_FULL,  // like FISHER_FULL, but with the matrix inverted
    COVARIANCE_LLT,   // like COVARIANCE_FULL, but stores the lower Cholesky factor of the matrix
};
 
}}}} // lsst::meas::multifit::sampling

%include "lsst/meas/multifit/sampling/RandomEngine.h"

%ignore lsst::meas::multifit::sampling::MixtureDistribution::draw;
%ignore lsst::meas::multifit::sampling::MixtureDistribution::update;

%std_nodefconst_type(lsst::meas::multifit::sampling::MixtureComponent);
%template(MixtureComponentList) std::vector<lsst::meas::multifit::sampling::MixtureComponent>;

%include "lsst/meas/multifit/sampling/MixtureDistribution.h"

%ignore lsst::meas::multifit::sampling::IterativeImportanceSampler::getSamples;

%extend lsst::meas::multifit::sampling::IterativeImportanceSampler {

    SampleIterator _getSamples() const {
        return SampleIterator(self->getSamples().begin(), self->getSamples().end());
    }

    PyObject * getLastSample() const {
        if (self->getSamples().empty()) {
            Py_RETURN_NONE;
        }
        return returnTable(self->getSamples().back());
    }

    %pythoncode {
        def getSamples(self):
            return self._getSamples()
    }
}

%include "lsst/meas/multifit/sampling/IterativeImportanceSampler.h"
