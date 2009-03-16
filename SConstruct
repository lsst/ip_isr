# -*- python -*-
#
# Setup our environment
#
import glob, os.path
import lsst.SConsUtils as scons

env = scons.makeEnv(
    "ip_isr",
    r"$HeadURL$",
    [
        ["boost", "boost/version.hpp", "boost_system:C++"],
        ["python", "Python.h"],
        ["cfitsio", "fitsio.h", "m cfitsio", "ffopen"],
        ["wcslib", "wcslib/wcs.h", "m wcs"], # remove m once SConsUtils bug fixed
        ["minuit", "Minuit/FCNBase.h", "lcg_Minuit:C++"],
        ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
        ["utils", "lsst/utils/Utils.h", "utils:C++"],
        ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
        ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
        ["security", "lsst/security/Security.h", "security:C++"],
        ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
        ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
        ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
        ["afw", "lsst/afw.h", "afw:C++"],
        ["meas_algorithms", "lsst/meas/algorithms/Interp.h", "meas_algorithms:C++"],
    ],
)
env.libs["ip_isr"] = env.getlibs("boost wcslib cfitsio minuit utils daf_base daf_data daf_persistence") \
    + env.getlibs("pex_exceptions pex_logging pex_policy security minuit afw meas_algorithms") \
    + env.libs["ip_isr"]

#
# Build/install things
#
for d in Split("doc include/lsst/ip/isr lib python/lsst/ip/isr tests examples"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [
    env.Install(env['prefix'], "python"),
    env.Install(env['prefix'], "include"),
    env.Install(env['prefix'], "lib"),
    env.InstallAs(os.path.join(env['prefix'], "doc", "doxygen"), os.path.join("doc", "htmlDir")),
    env.InstallEups(os.path.join(env['prefix'], "ups"), glob.glob(os.path.join("ups", "*.table"))),
    env.Install(env['prefix'], "pipeline"),
])

scons.CleanTree(r"*~ core *.so *.os *.o")

files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST Instrument Signature Removal package
""")

