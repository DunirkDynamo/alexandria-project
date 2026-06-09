"""
Analyzer-Plotter Wrapper Classes

These wrapper classes pair each analyzer with its corresponding plotter,
providing a convenient all-in-one interface for analysis and visualization.
"""

from .ctp401_wrapper import CTP401ModuleReporter
from .ctp404_wrapper import CTP404ModuleReporter
from .ctp515_wrapper import CTP515ModuleReporter
from .high_contrast_wrapper import HighContrastModuleReporter
from .uniformity_wrapper import UniformityModuleReporter

__all__ = [
    "UniformityModuleReporter",
    "HighContrastModuleReporter",
    "CTP401ModuleReporter",
    "CTP404ModuleReporter",
    "CTP515ModuleReporter",
]
