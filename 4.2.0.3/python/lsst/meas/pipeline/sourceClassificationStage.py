# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import itertools
import pdb
import sys

import lsst.pex.harness.stage as harnessStage
from lsst.pex.logging import Log, LogRec, endr
import lsst.pex.policy as pexPolicy
from lsst.meas.utils.sourceClassifier import SourceClassifier

__all__ = ["SourceClassificationStage", "SourceClassificationStageParallel"]

class SourceClassificationStageParallel(harnessStage.ParallelProcessing):
    """
    This stage implements source classification on pairs of sources.
    Two lists of sources, each corresponding to measurements from one
    of two per-visit exposures are expected to be passed in on the clipboard.
    These lists are expected to have identical size and order, such that the
    sources at position i in both lists are the pair of measurements
    from a single detection on the addition of a visits two difference images.
    The clipboard key names of these two lists is configurable via policy.

    The stage policy is additionally expected to contain a "Classifiers" key
    that describes up to 64 classifiers (which set source flag bits). Each of
    these must contain the following child keys:
    - "bits"         array of flag bits computed by the classifier (0-63)
    - "pythonClass"  the name of the python class for the classifier
    - "arguments"    "first":  process the first source list only
                     "second": process the second source list only
                     "both":   process both source lists
                     "pairs":  process pairs of sources from both lists
    - "parameters"   contains classifier specific configuration parameters 

    Note that classifiers are run in the order in which they appear in the
    stage policy file.

    See pipeline/SourceClassificationStageDictionary.paf for the details
    (including classifier specific parameters).

    Clipboard Input:
    - First difference source list (name determined by "sourceList1ClipboardKey" policy key)
    - Second difference source list (name determined by "sourceList2ClipboardKey" policy key)

    Clipboard Output:
    - the input clipboard is passed to the next stage with no structural
      modifications. Only individual sources are modified by this stage.
    """

    class Classifier(object):
        """
        Helper class for importing and running Ssource classifiers.
        """
        def __init__(self, policy, log):
            bits = policy.getIntArray("bits")
            pythonClass = policy.getString("pythonClass")
            components = pythonClass.split('.')
            className = components.pop()
            if len(components) == 0:
                raise RuntimeError("SourceClassifier package must be fully specified")
            else:
                try:
                    moduleName = '.'.join(components)
                    __import__(moduleName)
                    self._pythonClass = sys.modules[moduleName].__dict__[className]
                except Exception, e:
                    raise RuntimeError("Failed to instantiate class for SourceClassifier %s" % pythonClass)
            self._scope = policy.getString("scope")
            subpolicy = policy.getPolicy("parameter") if policy.exists("parameter") else Policy()
            if not issubclass(self._pythonClass, SourceClassifier):
                raise RuntimeError("%s is not a subclass of SourceClassifier - check stage policy" % pythonClass)
            # Create classifier, specifying a flag bit position
            self._classifier = self._pythonClass(bits, subpolicy)
            rec = LogRec(log, Log.INFO)
            rec << "Registered SourceClassifier" << { "pythonClass": pythonClass }
            for bit in bits:
                rec << { "bit": bit }
            rec << endr

        def invoke(self, sourceList1, sourceList2):
            if self._scope == "first":
                for source in sourceList1:
                    self._classifier.classify(source)
            elif self._scope == "second":
                for source in sourceList2:
                    self._classifier.classify(source)
            elif self._scope == "both":
                for source in itertools.chain(sourceList1, sourceList2):
                    self._classifier.classify(source)
            elif self._scope == "pairs":
                if sourceList1.size() != sourceList2.size():
                    raise RuntimeError("Source lists passed to classifier must have identical length")
                for s0, s1 in itertools.izip(sourceList1, sourceList2):
                    self._classifier.classify(s0, s1)
            else:
                raise RuntimeError("Unsupported source classifier argument type - check stage policy")

        def finish(self, log):
            self._classifier.finish(log)


    def setup(self):
        self.log = Log(self.log, "SourceClassificationStage - parallel")

        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
            "SourceClassificationStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        """
        Classify sources in the worker process
        """
        sourceList1 = clipboard.get(self.policy.getString("sourceList0Key"))
        sourceList2 = clipboard.get(self.policy.getString("sourceList1Key"))
        classifierList = []
        if self.policy.exists("classifier"):            
            classifierPolicyList = self.policy.getPolicyArray("classifier")      
            for classifierPolicy in classifierPolicyList:
                classifier = SourceClassificationStageParallel.Classifier(
                    classifierPolicy, self.log)
                classifierList.append(classifier)

        for classifier in classifierList:
            classifier.invoke(sourceList1, sourceList2)
        for classifier in classifierList:
            classifier.finish(self.log)

class SourceClassificationStage(harnessStage.Stage):
    parallelClass = SourceClassificationStageParallel

