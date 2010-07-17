# -*- python -*-
#
# Setup our environment
#
import glob, os.path, re, os
import lsst.SConsUtils as scons

try:
    scons.ConfigureDependentProducts
except AttributeError:
    import lsst.afw.SconsUtils
    scons.ConfigureDependentProducts = lsst.afw.SconsUtils.ConfigureDependentProducts

env = scons.makeEnv("meas_algorithms",
                    r"$HeadURL$",
                    scons.ConfigureDependentProducts("meas_algorithms"))

env.libs["meas_algorithms"] +=  env.getlibs("daf_base daf_data daf_persistence pex_logging pex_exceptions pex_policy security afw boost minuit2 utils wcslib")
#
# Build/install things
#
for d in Split("doc examples lib src tests") + \
        Split("include/lsst/meas/algorithms python/lsst/meas/algorithms"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", env.Install(env['prefix'], "examples"))
Alias("install", env.Install(env['prefix'], "include"))
Alias("install", env.Install(env['prefix'], "lib"))
Alias("install", env.Install(env['prefix'], "policy"))
Alias("install", env.Install(env['prefix'], "python"))
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

