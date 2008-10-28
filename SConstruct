# -*- python -*-
#
# Setup our environment
#
import glob, os.path
import lsst.SConsUtils as scons

env = scons.makeEnv(
    "ip_isr",
    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/ip/isr/trunk/SConstruct $",
    [
        ["boost", "boost/version.hpp", "boost_filesystem:C++"],
        ["vw", "vw/Math.h", "vwMath:C++"],
        ["python", "Python.h"],
        ["m", "math.h", "m", "sqrt"], # why is this optional?
        ["cfitsio", "fitsio.h", "m cfitsio", "ffopen"], # needed to link _isrLib.so
        ["wcslib", "wcslib/wcs.h", "m wcs"], # needed by afw
        ["minuit", "Minuit/FCNBase.h", "lcg_Minuit:C++"], # needed by afw
        ["lapack", None, "lapack", "dgesdd_"],
        ["utils", "lsst/utils/Utils.h", "utils:C++"],
        ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
        ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
        ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
        ["security", "lsst/security/Security.h", "security:C++"], # needed by daf_data
        ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
        ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"], # needed by daf_data
        ["daf_data", "lsst/daf/data.h", "daf_data:C++"], # needed by afw
        ["afw", "lsst/afw.h", "afw:C++"],
        ["detection", "lsst/detection/Footprint.h", "detection"],
    ],
)
env.libs["ip_isr"] = env.getlibs("boost vw lapack wcslib cfitsio utils daf_base pex_logging pex_exceptions pex_logging daf_persistence daf_data pex_policy minuit afw detection") + env.libs["ip_isr"]
env.libs["ip_isr"] += ["lapack"]

#
# Build/install things
#
for d in Split("doc include/lsst/ip/isr lib python/lsst/ip/isr tests examples"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", env.Install(env['prefix'], "python"))
Alias("install", env.Install(env['prefix'], "include"))
Alias("install", env.Install(env['prefix'], "lib"))
Alias("install", env.Install(env['prefix'], "pipeline"))
Alias("install", env.Install(env['prefix'] + "/bin", glob.glob("bin/*.py")))
Alias("install", env.InstallEups(env['prefix'] + "/ups", glob.glob("ups/*.table")))

scons.CleanTree(r"*~ core *.so *.os *.o")

files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST ISR packages
""")

