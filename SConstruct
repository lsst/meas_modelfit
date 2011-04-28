# -*- python -*-
#
# Setup our environment
#
# Do not change these
import glob, os.path, re, sys, traceback
import lsst.SConsUtils as scons


try:
    scons.ConfigureDependentProducts
except AttributeError:
    import lsst.afw.scons.SconsUtils
    scons.ConfigureDependentProducts = lsst.afw.scons.SconsUtils.ConfigureDependentProducts

env = scons.makeEnv(
    # The name of your package goes here.
    "meas_multifit",
    # This is used to try to get some version information.
    r"$HeadURL$",
    scons.ConfigureDependentProducts("meas_multifit")
)



# Describe what the package contains here.
env.Help("""
LSST Multifit Implementation.
""")

#
# Libraries needed to link libraries/executables
#
env.libs["meas_multifit"] += env.getlibs("boost wcslib cfitsio minuit2 gsl utils " +
        "daf_base daf_data daf_persistence pex_exceptions pex_logging pex_policy " +
        "security fftw3 afw ndarray")

if True:
    #
    # Workaround SConsUtils failure to find numpy .h files. Fixed in sconsUtils >= 3.3.2
    #
    import numpy
    env.Append(CCFLAGS = ["-I", numpy.get_include()])

# Build/install things
#
for d in (
    ".",
    "lib",
    "python/lsst/meas/multifit",
    "python/lsst/meas/multifit/sampling",
    #"examples",
    "tests",
    "doc",
):  
    if d != ".":
        try:
            SConscript(os.path.join(d, "SConscript"))
        except Exception, e:
            print >> sys.stderr, "In processing file %s:" % (os.path.join(d, "SConscript"))
            print >> sys.stderr, traceback.format_exc()
        Clean(d, Glob(os.path.join(d, "*~")))
        Clean(d, Glob(os.path.join(d, "*.pyc")))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [
    env.Install(env['prefix'], "include"),
    env.Install(env['prefix'], "lib"),
    env.Install(env['prefix'], "python"),
    env.Install(env['prefix'], "etc"),
    env.Install(env['prefix'], "policy"),
    env.Install(env['prefix'], "examples"),
    env.Install(env['prefix'], "tests"),
    env.Install(env['prefix'], "src"),
    env.InstallAs(os.path.join(env['prefix'], "doc", "doxygen"),
            os.path.join("doc", "htmlDir")),
    env.InstallEups(os.path.join(env['prefix'], "ups")),
])

scons.CleanTree(r"*~ core *.so *.os *.o")

#
# Build TAGS files
#
files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
