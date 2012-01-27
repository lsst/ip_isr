from lsst.ip.isr import IsrTask
from lsst.ip.isr import IsrTaskConfig
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

def makeCalibDict(butler, dataId):
    ret = {}
    for name in ("flat", "bias", "dark"):
        ret[name] = butler.get(name, dataId)
    return ret

dataid = {'visit':885449191, 'filter':'i', 'snap':0, 'raft':'2,2', 'sensor':'1,1', 'channel':'0,0'}
args = ['lsstSim', '/lsst3/weekly/data/obs_imSim-2011-09-07/PT1.2/', '--calib=/lsst3/weekly/data/obs_imSim-2011-09-07/PT1.2/', '--output=.']
parser = pipeBase.ArgumentParser()
namespace = parser.parse_args(argv=args, config=IsrTaskConfig())
calibSet = makeCalibDict(namespace.butler, dataid)
isrtask = IsrTask(config=namespace.config)
output = isrtask.run(namespace.butler.get("raw", dataid), calibSet)

