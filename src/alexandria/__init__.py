"""
Alexandria - Unified CatPhan Phantom Analysis Library

A comprehensive library for analyzing CatPhan CT phantom DICOM images.
Provides modular analyzers for uniformity, resolution, contrast, and linearity analysis.
"""

try:
    # `setuptools_scm` will write the current version to `_version.py` during build
    from ._version import __version__  # type: ignore
except Exception:
    # Fallback when package isn't installed via setuptools (dev/editable installs)
    __version__ = "0+unknown"

# Core analyzers (imported lazily by subpackages)
from .analyzers.uniformity import UniformityAnalyzer
from .analyzers.high_contrast import HighContrastAnalyzer
from .analyzers.ctp401 import CTP401Analyzer
from .analyzers.ctp404 import CTP404Analyzer
from .analyzers.ctp515 import CTP515Analyzer
from .analyzers.detailed_uniformity import DetailedUniformityAnalyzer

# Convenience wrappers (analyzer + plotter)
from .wrappers.uniformity_wrapper import UniformityModuleReporter
from .wrappers.high_contrast_wrapper import HighContrastModuleReporter
from .wrappers.ctp401_wrapper import CTP401ModuleReporter
from .wrappers.ctp404_wrapper import CTP404ModuleReporter
from .wrappers.ctp515_wrapper import CTP515ModuleReporter

__all__ = [
    'UniformityAnalyzer',
    'HighContrastAnalyzer',
    'CTP401Analyzer',
    'CTP404Analyzer',
    'CTP515Analyzer',
    'DetailedUniformityAnalyzer',
    'UniformityModuleReporter',
    'HighContrastModuleReporter',
    'CTP401ModuleReporter',
    'CTP404ModuleReporter',
    'CTP515ModuleReporter',
]
