from lsst.ip.isr import IsrTask
from lsst.ip.isr import IsrTaskConfig
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

args = ['lsstSim', '/lsst3/weekly/data/obs_imSim-2011-09-07/PT1.2/', '--calib=/lsst3/weekly/data/obs_imSim-2011-09-07/PT1.2/', '--output=/nfs/lsst/home/krughoff/lsst/lsst_devel/Linux64/LSST/DMS/ip/isr/examples/output/',\
         '--id', 'visit=885449191', 'snap=0', 'raft=2,2', 'sensor=1,1']
parser = pipeBase.ArgumentParser()
namespace = parser.parse_args(args=args, config=IsrTaskConfig())
isrtask = IsrTask(config=namespace.config)
print namespace.dataRefList
output = isrtask.run(namespace.dataRefList[0])

