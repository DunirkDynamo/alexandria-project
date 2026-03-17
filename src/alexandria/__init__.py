"""
Alexandria - Unified CatPhan Phantom Analysis Library

A comprehensive library for analyzing CatPhan CT phantom DICOM images.
Provides modular analyzers for uniformity, resolution, contrast, and linearity analysis.
"""

__version__ = "1.0.0"

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
