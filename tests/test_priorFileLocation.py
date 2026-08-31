#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import os
import unittest

import lsst.utils.tests
import lsst.meas.modelfit


class PriorFileLocationTestCase(lsst.utils.tests.TestCase):
    """Test that FILE priors are found both when the MEAS_MODELFIT_DIR
    environment variable is set (EUPS-provided) and when it is not set, in
    which case the data directory is located relative to the meas_modelfit
    shared library itself.
    """

    def setUp(self):
        # Remember the current value so we can restore it in tearDown.
        self.savedDir = os.environ.get("MEAS_MODELFIT_DIR")
        # Determine a valid package directory to point the environment
        # variable at for the "env var set" case. Prefer the value already in
        # the environment; otherwise ask EUPS/utils for it.
        if self.savedDir is not None:
            self.packageDir = self.savedDir
        else:
            try:
                self.packageDir = lsst.utils.getPackageDir("meas_modelfit")
            except LookupError:
                self.packageDir = None

    def tearDown(self):
        if self.savedDir is not None:
            os.environ["MEAS_MODELFIT_DIR"] = self.savedDir
        else:
            os.environ.pop("MEAS_MODELFIT_DIR", None)

    def _makeControl(self):
        ctrl = lsst.meas.modelfit.CModelStageControl()
        ctrl.priorSource = "FILE"
        ctrl.priorName = "s13-v2-disk-08"
        return ctrl

    def testFindsPriorWithEnvVar(self):
        """The prior is found via the MEAS_MODELFIT_DIR environment variable."""
        if self.packageDir is None:
            self.skipTest("Cannot determine meas_modelfit package directory")
        os.environ["MEAS_MODELFIT_DIR"] = self.packageDir
        prior = self._makeControl().getPrior()
        self.assertIsNotNone(prior)

    def testFindsPriorWithoutEnvVar(self):
        """The prior is found relative to the shared library when the
        environment variable is not set.
        """
        os.environ.pop("MEAS_MODELFIT_DIR", None)
        prior = self._makeControl().getPrior()
        self.assertIsNotNone(prior)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
