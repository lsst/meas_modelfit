// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "ndarray/eigen.h"

#include "lsst/pex/exceptions.h"
#define LSST_MAX_DEBUG 0
#include "lsst/pex/logging/Debug.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/meas/modelfit/psf.h"
#define LSST_DEBUG_DUMP
#include "lsst/meas/modelfit/DebugDump.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

// Create a column-major matrix (ndarray creates row-major matrices directly)
ndarray::Array<Pixel,2,-2> makeMatrix(int nRows, int nCols) {
    return ndarray::Array<Pixel,2,2>(ndarray::allocate(nCols, nRows)).transpose();
}

} // anonymous

class PsfFitter::Impl {
public:

    Impl(afw::geom::Box2I const & bbox, std::vector<PsfFitterComponentControl> const & ctrls) :
        bbox(bbox),
        _x(ndarray::allocate(bbox.getArea())),
        _y(ndarray::allocate(bbox.getArea()))
    {
        int const width = bbox.getWidth();
        int const height = bbox.getHeight();

        // Fill flattened 1-d arrays with positions
        auto xIter = _x.begin();
        auto yIter = _y.begin();
        for (int i = 0, y = bbox.getMinY(); i < height; ++i, ++y) {
            for (int j = 0, x = bbox.getMinX(); j < width; ++j, ++x) {
                *xIter = x;
                *yIter = y;
                ++xIter;
                ++yIter;
            }
        }

        // Create factories for builders that evaluate Gaussian and Shapelet functions.
        // We construct all the factories first so we can compute the workspace needed
        // by the builders, so the builders can share that workspace.
        // See shapelet::MatrixBuilderFactory for more info.
        std::vector<shapelet::MatrixBuilderFactory<Pixel>> factories;
        factories.reserve(ctrls.size() + 1);
        factories.push_back(shapelet::MatrixBuilderFactory<Pixel>(_x, _y, 0));
        int workspaceSize = factories.back().computeWorkspace();
        for (auto const & ctrl : ctrls) {
            factories.push_back(shapelet::MatrixBuilderFactory<Pixel>(_x, _y, ctrl.order));
            workspaceSize = std::max(workspaceSize, factories.back().computeWorkspace());
        }
        // Construct the shared workspace.
        shapelet::MatrixBuilderWorkspace<Pixel> workspace(workspaceSize);
        // Construct the builders using the shared workspace.
        _builders.reserve(factories.size());
        for (auto const & factory : factories) {
            shapelet::MatrixBuilderWorkspace<Pixel> wsCopy(workspaceSize);
            _builders.push_back(factory(wsCopy));
        }
    }

    void evaluateGaussian(
        ndarray::Array<Pixel,2,-1> const & matrix,
        afw::geom::ellipses::Ellipse const & ellipse
    ) const {
        matrix.deep() = 0.0;
        return _builders.front()(matrix, ellipse);
    }

    void evaluateShapelet(
        ndarray::Array<Pixel,2,-1> const & matrix,
        afw::geom::ellipses::Ellipse const & ellipse,
        int index
    ) const {
        matrix.deep() = 0.0;
        return _builders[index + 1](matrix, ellipse);
    }

    // Create an ellipse core, after fudging it to make sure it isn't singular.
    afw::geom::ellipses::Quadrupole makeSafeEllipse(double ixx, double iyy, double ixy, double rMin) const {
        // We start by "deconvolving" a circle with the min radius (subtracting it in quadrature).
        ixx -= rMin*rMin;
        iyy -= rMin*rMin;
        // Now we "clip" the ellipse - set ixx and iyy and the determinant to zero if negative.
        ixx = std::max(ixx, 0.0);
        iyy = std::max(ixx, 0.0);
        if (ixx*iyy < ixy*ixy) {
            ixy = std::sqrt(ixx*iyy);
        }
        // We now "convolve" by a circle with the min radius.  The result is a no-op if the
        // ellipse was already strictly larger than the min radius circle, and a well-behaved
        // similar ellipse if it wasn't.
        ixx += rMin*rMin;
        iyy += rMin*rMin;
        return afw::geom::ellipses::Quadrupole(ixx, iyy, ixy);
    }

    int fitEM(
        shapelet::MultiShapeletFunction & result,
        ndarray::Array<Pixel const ,1,1> const & data,
        std::vector<PsfFitterComponentControl> const & ctrls,
        int maxIterations,
        double absTol,
        double relTol
    ) const {

        pex::logging::Debug log("meas.modelfit.PsfFitter");
        DebugDump dump("psf.fitEM");

        std::size_t const nComponents = result.getComponents().size();
        std::size_t const nPix = data.getSize<0>();

        ndarray::Array<Pixel,2,-2> matrix = makeMatrix(nPix, nComponents);
        Eigen::Matrix<Pixel,Eigen::Dynamic,1> model(nPix);
        Eigen::Array<Pixel,Eigen::Dynamic,1> weighted(nPix);
        Eigen::Array<Pixel,Eigen::Dynamic,1> dx(nPix);
        Eigen::Array<Pixel,Eigen::Dynamic,1> dy(nPix);

        // Best-fit coefficients of the linear fit from last iteration.
        // Initialize with the 0th-order coefficients from the input MultiShapeletFunction.
        ndarray::Array<double,1,1> lastCoeff = ndarray::allocate(nComponents);
        for (std::size_t n = 0; n < nComponents; ++n) {
            lastCoeff[n] = result.getComponents()[n].getCoefficients()[0];

            log.debug<7>(
                "Initial state for component %d: coefficient=%g, moments=(%g, %g, %g), center=(%g, %g)",
                n, lastCoeff[n],
                afw::geom::ellipses::Quadrupole(result.getComponents()[n].getEllipse().getCore()).getIxx(),
                afw::geom::ellipses::Quadrupole(result.getComponents()[n].getEllipse().getCore()).getIyy(),
                afw::geom::ellipses::Quadrupole(result.getComponents()[n].getEllipse().getCore()).getIxy(),
                result.getComponents()[n].getEllipse().getCenter().getX(),
                result.getComponents()[n].getEllipse().getCenter().getY()
            );

        }

        // Initialize a nonlinear least squares fitter so we can reuse it within the loop.
        afw::math::LeastSquares lstsq(afw::math::LeastSquares::NORMAL_EIGENSYSTEM, nComponents);

        // Expectation-Maximization loop.  We update the ellipses of the result MultiShapeletFunction in-place,
        // but don't do anything with the coefficients yet, as whatever we might have done would have been
        // totally overwritten by the final linear fit.
        for (int iteration = 0; iteration < maxIterations; ++iteration) {

            // ----------------------------------------------------------------------------------------------
            // Expectation
            // Do a linear fit (holding ellipses fixed) of a set of Gaussians (one for each shapelet
            // component) to the image.
            // ----------------------------------------------------------------------------------------------

            // Evaluate the model matrix, column by column.  Each column is a Gaussian corresponding
            // to one of the component ellipses, normalized to the ShapeletFunction convention.
            for (std::size_t n = 0; n < nComponents; ++n) {
                evaluateGaussian(matrix[ndarray::view()(n,n+1)], result.getComponents()[n].getEllipse());
                log.debug<7>(
                    "Iteration %d for component %d: matrix column norm is %g",
                    iteration, n, matrix[ndarray::view()(n)].asEigen().norm()
                );
            }
            dump.copy((boost::format("iter=%05d:matrix") % iteration).str(), matrix);

            // Solve for the best-fit linear combination of components.
            lstsq.setDesignMatrix(matrix, data);

            dump.copy((boost::format("iter=%05d:solution") % iteration).str(), lstsq.getSolution());

            // Compute the model corresponding to the best fit solution.
            model = matrix.asEigen() * lstsq.getSolution().asEigen().cast<Pixel>();  // vector = matrix*vector

            dump((boost::format("iter=%05d:model") % iteration).str(), model);

            // ----------------------------------------------------------------------------------------------
            // Maximization
            // ----------------------------------------------------------------------------------------------

            for (std::size_t n = 0; n < nComponents; ++n) {
                // Weight the data image by relative contribution of each component to each pixel in the
                // best-fit model
                if (lstsq.getSolution()[n] < 0) {
                    // Negative components aren't allowed to contribute; we just don't update the ellipse
                    continue;
                }
                weighted =
                    data.asEigen<Eigen::ArrayXpr>()
                    * lstsq.getSolution()[n]
                    * matrix.asEigen<Eigen::ArrayXpr>().col(n)
                    / model.array();
                // Compute moments and update the ellipse
                dx = _x.asEigen<Eigen::ArrayXpr>();
                dy = _y.asEigen<Eigen::ArrayXpr>();
                double i0 = weighted.sum();
                double ix = (weighted * dx).sum() / i0;
                double iy = (weighted * dy).sum() / i0;
                result.getComponents()[n].getEllipse().setCenter(afw::geom::Point2D(ix, iy));
                dx -= ix;
                dy -= iy;
                double ixx = (weighted * dx.square()).sum() / i0;
                double iyy = (weighted * dy.square()).sum() / i0;
                double ixy = (weighted * dx * dy).sum() / i0;
                result.getComponents()[n].getEllipse().setCore(
                    makeSafeEllipse(ixx, iyy, ixy, ctrls[n].radiusMin)
                );
                log.debug<7>(
                    "Iteration %d for component %d: coefficient=%g, moments=(%g, %g, %g), center=(%g, %g)",
                    iteration, n, lstsq.getSolution()[n], ixx, iyy, ixy, ix, iy
                );
            }

            // Check whether the linear fit results changed significantly in the last iteration; if not,
            // declare convergence and quit.
            bool converged = true;
            for (std::size_t n = 0; n < nComponents; ++n) {
                double absDiff = std::abs(lstsq.getSolution()[n] - lastCoeff[n]);
                double absMean = 0.5*(std::abs(lstsq.getSolution()[n]) + std::abs(lastCoeff[n]));
                if (absDiff > absTol || absDiff > absMean*relTol) {
                    converged = false;
                    break;
                }
            }
            if (converged) return iteration;
        }
        return maxIterations;
    }

    void fitShapelets(
        shapelet::MultiShapeletFunction & result,
        ndarray::Array<Pixel const ,1,1> const & data
    ) const {
        std::size_t const nComponents = result.getComponents().size();
        std::size_t const nPix = data.getSize<0>();

        // Figure out how big the matrix will need to be.
        std::size_t nTotalCoefficients = 0;
        for (std::size_t n = 0; n < nComponents; ++n) {
            nTotalCoefficients += result.getComponents()[n].getCoefficients().size();
        }
        ndarray::Array<Pixel,2,-2> matrix = makeMatrix(nPix, nTotalCoefficients);
        // Evaluate the model matrix, in chunks of columns corresponding to a single shapelet component.
        std::size_t offset = 0;
        for (std::size_t n = 0; n < nComponents; ++n) {
            int nCols = result.getComponents()[n].getCoefficients().getSize<0>();
            evaluateShapelet(matrix[ndarray::view()(offset, offset + nCols)],
                             result.getComponents()[n].getEllipse(), n);
            offset += nCols;
        }
        // Do the actual fit.
        auto lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, data);
        // Put the fit results back into the MultiShapeletFunction.
        offset = 0;
        for (std::size_t n = 0; n < nComponents; ++n) {
            int nCols = result.getComponents()[n].getCoefficients().getSize<0>();
            result.getComponents()[n].getCoefficients().deep()
                = lstsq.getSolution()[ndarray::view(offset, offset+nCols)];
            offset += nCols;
        }
    }

    afw::geom::Box2I const bbox;

private:
    ndarray::Array<Pixel,1,1> _x;
    ndarray::Array<Pixel,1,1> _y;
    std::vector<shapelet::MatrixBuilder<Pixel>> _builders;
};

PsfFitter::PsfFitter(PsfFitterControl const & ctrl) :
    _maxIterations(ctrl.maxIterations),
    _absTol(ctrl.absTol),
    _relTol(ctrl.relTol),
    _innerIndex(-1),
    _primaryIndex(-1),
    _wingsIndex(-1),
    _outerIndex(-1)
{
    int n = 0;
    if (ctrl.inner.order >= 0) {
        _ctrls.push_back(ctrl.inner);
        _innerIndex = n;
        ++n;
    }
    if (ctrl.primary.order >= 0) {
        _ctrls.push_back(ctrl.primary);
        _primaryIndex = n;
        ++n;
    } else {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "PsfFitter control must have a primary component with nonnegative order"
        );
    }
    if (ctrl.wings.order >= 0) {
        _ctrls.push_back(ctrl.wings);
        _wingsIndex = n;
        ++n;
    }
    if (ctrl.outer.order >= 0) {
        _ctrls.push_back(ctrl.outer);
        _outerIndex = n;
    }
    // _impl is constructed on first use
}

shapelet::MultiShapeletFunctionKey PsfFitter::addModelFields(
    afw::table::Schema & schema,
    std::string const & prefix
) const {
    std::vector<int> orders;
    orders.reserve(_ctrls.size());
    for (auto const & ctrl : _ctrls) {
        orders.push_back(ctrl.order);
    }
    return shapelet::MultiShapeletFunctionKey::addFields(
        schema, prefix, "multi-Shapelet approximation to the PSF model",
        "pixels", // ellipse units
        "",       // coefficient units (unitless)
        orders
    );
}


void PsfFitter::adaptComponent(
    shapelet::MultiShapeletFunction & current,
    shapelet::MultiShapeletFunction const & previous,
    PsfFitter const & previousFitter,
    int currentIndex,
    int previousIndex
) const {
    if (currentIndex < 0) return;
    shapelet::ShapeletFunction currentComponent(_ctrls[currentIndex].order, shapelet::HERMITE);
    if (previousIndex < 0) {
        auto const & previousPrimary = previous.getComponents()[previousFitter._primaryIndex];
        // previous didn't include the component we're setting up now, so we base its parameters
        // off the previous primary component (but we initialize the coefficients to zero).
        currentComponent.setEllipse(previousPrimary.getEllipse());
        currentComponent.getEllipse().getCore().scale(_ctrls[currentIndex].radiusFactor);
    } else {
        // previous did include the component we're setting up now, so we copy it over,
        // including as many coeffients as possible.
        auto const & previousComponent = previous.getComponents()[previousIndex];
        currentComponent.setEllipse(previousComponent.getEllipse());
        int minOrder = std::min(previousFitter._ctrls[previousIndex].order, _ctrls[currentIndex].order);
        int minSize = shapelet::computeSize(minOrder);
        currentComponent.getCoefficients()[ndarray::view(0, minSize)]
            = previousComponent.getCoefficients()[ndarray::view(0, minSize)];
    }
    current.getComponents().push_back(std::move(currentComponent));
}

shapelet::MultiShapeletFunction PsfFitter::adapt(
    shapelet::MultiShapeletFunction const & previous,
    PsfFitter const & previousFitter
) const {
    shapelet::MultiShapeletFunction current;
    adaptComponent(current, previous, previousFitter, _innerIndex, previousFitter._innerIndex);
    adaptComponent(current, previous, previousFitter, _primaryIndex, previousFitter._primaryIndex);
    adaptComponent(current, previous, previousFitter, _wingsIndex, previousFitter._wingsIndex);
    adaptComponent(current, previous, previousFitter, _outerIndex, previousFitter._outerIndex);
    return current;
}


shapelet::MultiShapeletFunction PsfFitter::makeInitial(
    afw::geom::ellipses::Quadrupole const & moments
) const {
    shapelet::MultiShapeletFunction result;
    for (auto const & ctrl : _ctrls) {
        result.getComponents().push_back(shapelet::ShapeletFunction(ctrl.order, shapelet::HERMITE));
        shapelet::ShapeletFunction & current = result.getComponents().back();
        current.getEllipse().setCore(moments);
        current.getEllipse().getCore().scale(ctrl.radiusFactor);
    }
    result.getComponents()[_primaryIndex].getCoefficients()[0] = 1.0 / shapelet::ShapeletFunction::FLUX_FACTOR;
    return result;
}

namespace {

// Create a 1-d array from an image by flattening over columns and then rows.
ndarray::Array<Pixel const,1,1> flattenImage(afw::image::Image<Pixel> const & image) {
    // If the image is already contiguous, we can just make a view with different shape and strides.
    ndarray::Array<Pixel const,2,2> contiguous = ndarray::dynamic_dimension_cast<2>(image.getArray());
    if (contiguous.isEmpty()) { // If it isn't contiguous, we have to deep-copy.
        contiguous = ndarray::copy(image.getArray());
    }
    return ndarray::flatten<1>(contiguous);
}

} // anonymous

int PsfFitter::apply(
    shapelet::MultiShapeletFunction & model,
    afw::image::Image<Pixel> const & image
) const {

    // Construct an Impl object for cached PSF-bbox-dependent state if necessary.
    if (!_impl || _impl->bbox != image.getBBox()) {
        _impl.reset(new Impl(image.getBBox(), _ctrls));
    }

    // View or copy of the image in 1-d, shape=(height*width)
    ndarray::Array<Pixel const,1,1> data1d = flattenImage(image);

    // Expectation-Maximization loop.  We update the ellipses of the result MultiShapeletFunction in-place,
    // but don't do anything with the coefficients yet, as whatever we might have done would have been
    // totally overwritten by the final linear fit.
    int nIterations = _impl->fitEM(model, data1d, _ctrls, _maxIterations, _absTol, _relTol);

    // Do a final linear fit with all shapelet terms active.
    _impl->fitShapelets(model, data1d);

    return nIterations;
}

PsfFitter::~PsfFitter() {} // needs to be in .cc for forward-declared Impl

}}} // namespace lsst::meas::modelfit
