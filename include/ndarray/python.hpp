#ifndef NDARRAY_python_hpp_INCLUDED
#define NDARRAY_python_hpp_INCLUDED

/**
 *  @file ndarray/python.hpp
 *  \brief Public header file for ndarray Python support.
 *
 *  \warning Both the main Python C-API header, "Python.h", and the
 *  Numpy C-API headers "arrayobject.h" and "ufuncobject.h" must
 *  be included before ndarray/python.hpp or any of the files in
 *  ndarray/python.
 *
 *  \note This file is not included by the main "ndarray.hpp" header file.
 */

/** \defgroup PythonGroup Python Support
 *
 *  The ndarray Python support module provides conversion
 *  functions between ndarray objects, notably Array and
 *  Vector, and Python Numpy objects.
 *
 *  \note The Numpy C-API header files must be included
 *  <em>before</em> ndarray-python.hpp, and the
 *  <b><tt>import_array()</tt></b> and
 *  <b><tt>import_ufunc()</tt></b> functions must be
 *  called appropriately (see the Numpy C-API documentation
 *  for more details).
 */

/// \defgroup PythonInternalGroup Python Support Internals

/// \internal @namespace ndarray::detail \brief Internal namespace for ndarray Python support.

#include "ndarray.hpp"
#include "ndarray/python/numpy.hpp"
#include "ndarray/python/ufunctors.hpp"
#include "ndarray/python/Vector.hpp"

#endif // !NDARRAY_python_hpp_INCLUDED
