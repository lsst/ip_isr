from __future__ import with_statement
import sys
import os
import subprocess
import time
import lsst.pex.logging

def startPipeline(nodeList, pipelinePolicy, runId):
    """Start pipeline execution
    
    Inputs:
    - nodeList: path to mpi machine file; environment variables and relative paths are expanded
    
    The node list file must be in the pipeline directory; this directory also contains
    the "policy" directory (containing pipeline policy files) and must be writable.
    
    The node list file contains information about the nodes on which to run the pipeline.
    In its simplest form there is one line per node in the format:
       ipaddress[:nslices]
    where nslices is the number of CPUs that you wish to use on that node; the defaults is 1.
    Blank lines and lines beginning with # are ignored.
    Additional options may be specified; see documentation for MPICH2 used with
    the MPD process management environment.

    The pipeline uses one slice just to run the preprocess and postprocess phase;
    all other slices are used to run slices of the process phase.
    For example:
    - To run one slice of the process phase on one host (that has at least 2 CPUs):
      specify one host with 2 slices (1 for pre/postprocess and 1 for process)
    - To run three slices of the process phase on two hosts (each with at least 2 CPUs):
      specify two hosts with 2 slices each (4 slices: 1 for pre/postprocess and 3 for process)
      
    Note for those coming from run.sh: nslices in the code below includes the slice for
    the pre/postprocess phase and thus is one greater than nslices in run.sh;
    as a result the universe size = nslices.
    """
    nodeList = os.path.abspath(os.path.expandvars(nodeList))
    pipelineDir = os.path.dirname(nodeList)
    lsst.pex.logging.Trace("pex.harness.startPipeline", 3, "pipelineDir=%s" % (pipelineDir,))

    nnodes, nslices = parseNodeList(nodeList)
    lsst.pex.logging.Trace("pex.harness.startPipeline", 3, "nnodes=%s; nslices=%s" % (nnodes, nslices))
    
    lsst.pex.logging.Trace("pex.harness.startPipeline", 3, "Running mpdboot")
    subprocess.call(["mpdboot", "--totalnum=%s" % (nnodes,), "--file=%s" % (nodeList,), "--verbose"])
    
    time.sleep(3)
    lsst.pex.logging.Trace("pex.harness.startPipeline", 3, "Running mpdtrace")
    subprocess.call(["mpdtrace", "-l"])
    time.sleep(2)
    
    lsst.pex.logging.Trace("pex.harness.startPipeline", 3, "Running mpiexec")
    subprocess.call(
        ["mpiexec", "-usize", str(nslices), "-machinefile", nodeList, "-np", "1",
        "runPipeline.py", pipelinePolicy, runId],
        cwd = pipelineDir,
    )
    
    time.sleep(1)    
    lsst.pex.logging.Trace("pex.harness.startPipeline", 3, "Running mpdallexit")
    subprocess.call(["mpdallexit"])

def parseNodeList(nodeList):
    """Return (nnodes, nslices)
    where:
    - nnodes = number of nodes (host IP addresses) in file
    - nslices = total number of slices specified    
    """
    nnodes = 0
    nslices = 0
    with file(nodeList, "r") as nodeFile:
        for line in nodeFile:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue

            try:
                # strip optional extra arguments
                hostInfo = line.split()[0]
                # parse optional number of slices
                hostSlice = hostInfo.split(":")
                if len(hostSlice) == 1:
                    nslices += 1
                elif len(hostSlice) == 2:
                    nslices += int(hostSlice[1])
                else:
                    raise RuntimeError("Could not parse host info %r" % hostInfo)
            except Exception, e:
                raise RuntimeError("Cannot parse nodeList line %r; error = %s" % (line, e))
            nnodes += 1
    return (nnodes, nslices)
            
            
