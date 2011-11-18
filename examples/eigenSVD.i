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


%feature("autodoc", "1");
%module(package="eigenSVD") eigenSVD

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored
#pragma SWIG nowarn=401                 // nothin known about base class X
%{
#include "lsst/afw/detection.h"
#include "lsst/ndarray/eigen.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#define PY_ARRAY_UNIQUE_SYMBOL EIGEN_SVD_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "lsst/ndarray/python.h"
#include "lsst/ndarray/python/eigen.h"
%}

/******************************************************************************/
%init %{
    import_array();
%}

%include "lsst/ndarray/ndarray.i"

/*****************************************************************************/
%declareNumPyConverters(lsst::ndarray::Array<double,2,2>);
%declareNumPyConverters(lsst::ndarray::Array<double,1,1>);

%inline %{
    void _svd(
        lsst::ndarray::Array<double,2,2> const & m, 
        lsst::ndarray::Array<double,2,2> const & u, 
        lsst::ndarray::Array<double,1,1> const & s, 
        lsst::ndarray::Array<double,2,2> const & v
    ) {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd;
        svd.compute(m.asEigen(), Eigen::ComputeThinU | Eigen::ComputeThinV);
        u.asEigen() = svd.matrixU();
        s.asEigen() = svd.singularValues();
        v.asEigen() = svd.matrixV();
    }
    void _reconstruct(
        lsst::ndarray::Array<double,2,2> const & m, 
        lsst::ndarray::Array<double,2,2> const & u, 
        lsst::ndarray::Array<double,1,1> const & s, 
        lsst::ndarray::Array<double,2,2> const & v
    ) {
        m.asEigen() = u.asEigen() * s.asEigen().asDiagonal()
	    * v.asEigen().transpose();
    }
%}
    
%pythoncode %{
    def svd(m):
        import numpy
        u = numpy.zeros(m.shape, dtype=float)
        s = numpy.zeros(m.shape[1], dtype=float)
        v = numpy.zeros((m.shape[1], m.shape[1]), dtype=float)
        _svd(m, u, s, v)
        return u, s, v.transpose()
    def reconstruct(u, s, vt):
        import numpy
        m = numpy.zeros(u.shape, dtype=float)
        _reconstruct(m, u.copy(), s.copy(), vt.transpose().copy())
        return m
%}
