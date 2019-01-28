# This file is part of ip_isr.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys

from lsst.ip.isr import AssembleCcdTask
import lsst.afw.display as afwDisplay
import exampleUtils


def runAssembler():
    """Run the task to assemble amps into a ccd.
    """

    # Create the assemble task with default config
    assembleConfig = AssembleCcdTask.ConfigClass()
    assembleTask = AssembleCcdTask(config=assembleConfig)
    frame = 0

    # The assemble task can operate in two ways:
    # 1. On a dictionary of amp size images
    # 2. On a single amp mosaic image
    # Both methods should produce the same output.
    for isPerAmp in (True, False):
        assemblyInput = exampleUtils.makeAssemblyInput(isPerAmp)
        assembledExposure = assembleTask.assembleCcd(assemblyInput)
        afwDisplay.Display(frame=frame).mtv(assembledExposure.getMaskedImage(),
                                            title="Per amp input is %s" % isPerAmp)
        frame += 1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of AssembleCcdTask")
    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    args = parser.parse_args()

    if args.debug:
        try:
            import debug  # noqa F401
        except ImportError as e:
            print(e, file=sys.stderr)
    runAssembler()
