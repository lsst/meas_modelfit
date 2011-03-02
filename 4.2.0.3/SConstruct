# -*- python -*-
#
# Setup our environment
#
import os.path, re, os
import lsst.SConsUtils as scons

env = scons.makeEnv("meas_pipeline",
                    r"$HeadURL$",
                    []
                    )
#
# Build/install things
#
for d in Split("doc tests python/lsst/meas/pipeline"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", env.Install(env['prefix'], "python"))
Alias("install", env.Install(env['prefix'], "policy"))
Alias("install", env.InstallEups(os.path.join(env['prefix'], "ups")))

scons.CleanTree(r"*~ core *.so *.os *.o")
#
# Build TAGS files
#
files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST FrameWork packages
""")

