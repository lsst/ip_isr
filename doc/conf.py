"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.ip.isr


_g = globals()
_g.update(build_package_configs(
    project_name='ip_isr',
    version=lsst.ip.isr.version.__version__))
