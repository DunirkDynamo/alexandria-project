"""
High Contrast Module Wrapper

Combines HighContrastAnalyzer with HighContrastPlotter for convenient
MTF/resolution analysis and visualization.
"""

import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict, Any
from ..analyzers.high_contrast import HighContrastAnalyzer


class HighContrastModuleReporter:
    def __init__(
        self,
        image: Optional[Any] = None,
        pixel_spacing: Optional[float] = None,
        center: Optional[Tuple[float, float]] = None,
        t_offset_deg: float = 0.0,
        rotation_offset: Optional[float] = None,
        dicom_set: Optional[List] = None,
        slice_index: Optional[int] = None,
        lp_r_mm: float = 48.0,
        samples_per_segment: int = 50,
    ):
        self.analyzer = HighContrastAnalyzer(
            image=image,
            pixel_spacing=pixel_spacing,
            center=center,
            t_offset_deg=t_offset_deg,
            rotation_offset=rotation_offset,
            dicom_set=dicom_set,
            slice_index=slice_index,
            lp_r_mm=lp_r_mm,
            samples_per_segment=samples_per_segment,
        )
        self.results = None
        self.figure = None

    def analyze(self, write_log: bool = False, verbose: bool = True) -> Dict[str, Any]:
        self.results = self.analyzer.analyze(write_log=write_log, verbose=verbose)
        return self.results

    def plot(self, show: bool = False, **kwargs) -> plt.Figure:
        if self.results is None:
            self.analyze(verbose=False)
        from ..plotters.high_contrast_plotter import HighContrastPlotter
        plotter = HighContrastPlotter(self.analyzer)
        self.figure = plotter.plot(**kwargs)
        if show:
            plt.show()
        return self.figure

    def analyze_and_plot(self, verbose: bool = True, show: bool = False, **kwargs) -> Tuple[Dict[str, Any], plt.Figure]:
        results = self.analyze(verbose=verbose)
        figure = self.plot(show=show, **kwargs)
        return results, figure

    def save_plot(self, filepath: str, dpi: int = 300, **kwargs):
        if self.figure is None:
            self.plot()
        self.figure.savefig(filepath, dpi=dpi, bbox_inches='tight', **kwargs)

    def get_summary(self) -> Dict[str, str]:
        if self.results is None:
            self.analyze(verbose=False)
        return self.analyzer.get_results_summary()

    def close_plot(self):
        if self.figure is not None:
            plt.close(self.figure)
            self.figure = None
