/*
 * This file is part of meas_modelfit.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;
using lsst::cpputils::python::WrapperCollection;

namespace lsst {
namespace meas {
namespace modelfit {

void wrapAdaptiveImportanceSampler(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapCmodel(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapIntegrals(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapLikelihood(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapMixture(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapModel(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapMultiModel(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapOptimizer(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapPixelFitRegion(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapPriors(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapPsf(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapSampler(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapTruncatedGaussian(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapUnitSystem(lsst::cpputils::python::WrapperCollection &wrappers);
void wrapUnitTransformedLikelihood(lsst::cpputils::python::WrapperCollection &wrappers);

PYBIND11_MODULE(_modelfitLib, mod) {
    lsst::utils::python::WrapperCollection wrappers(mod, "lsst.meas.modelfit");

    wrappers.addInheritanceDependency("lsst.meas.base");

    wrappers.addSignatureDependency("lsst.afw.detection");
    wrappers.addSignatureDependency("lsst.afw.image");
    wrappers.addSignatureDependency("lsst.afw.math");
    wrappers.addSignatureDependency("lsst.afw.geom.ellipses");
    wrappers.addSignatureDependency("lsst.afw.table");
    wrappers.addSignatureDependency("lsst.shapelet");

    wrapPriors(wrappers);
    wrapUnitSystem(wrappers);
    wrapModel(wrappers);
    wrapLikelihood(wrappers);
    wrapMixture(wrappers);
    wrapOptimizer(wrappers);
    wrapSampler(wrappers);
    wrapPixelFitRegion(wrappers);

    wrapAdaptiveImportanceSampler(wrappers);
    wrapCmodel(wrappers);
    wrapIntegrals(wrappers);
    wrapMultiModel(wrappers);
    wrapPsf(wrappers);
    wrapTruncatedGaussian(wrappers);
    wrapUnitTransformedLikelihood(wrappers);
    wrappers.finish();
}

}
}
}
