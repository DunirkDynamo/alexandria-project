"""
CTP401 Module Wrapper

Combines CTP401Analyzer with CTP401Plotter for convenient
4-ROI linearity analysis and visualization in a single class.
"""

import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict, Any
from ..analyzers.ctp401 import CTP401Analyzer


class CTP401ModuleReporter:
    def __init__(
        self,
        image: Optional[Any] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,
        dicom_set: Optional[List[Any]] = None,
        slice_index: Optional[int] = None,
        roi_radius: float = 3.5,
        material_distance: float = 58.5,
        edge_threshold: float = 100.0,
    ):
        self.analyzer = CTP401Analyzer(
            image=image,
            center=center,
            pixel_spacing=pixel_spacing,
            spacing=spacing,
            dicom_set=dicom_set,
            slice_index=slice_index,
            roi_radius=roi_radius,
            material_distance=material_distance,
            edge_threshold=edge_threshold,
        )
        self.results = None
        self.figure = None

    def analyze(self, t_offset: float = 0.0, verbose: bool = True) -> Dict[str, Any]:
        self.results = self.analyzer.analyze(t_offset=t_offset, verbose=verbose)
        return self.results

    def detect_rotation(self, initial_angle_deg: float = 0.0) -> float:
        return self.analyzer.detect_rotation(initial_angle_deg=initial_angle_deg)

    def plot(self, show: bool = False, **kwargs) -> plt.Figure:
        if self.results is None:
            self.analyze(verbose=False)
        from ..plotters.ctp401_plotter import CTP401Plotter
        plotter = CTP401Plotter(self.analyzer)
        self.figure = plotter.plot(**kwargs)
        if show:
            plt.show()
        return self.figure

    def analyze_and_plot(self, t_offset: float = 0.0, verbose: bool = True, show: bool = False, **kwargs) -> Tuple[Dict[str, Any], plt.Figure]:
        results = self.analyze(t_offset=t_offset, verbose=verbose)
        fig = self.plot(show=show, **kwargs)
        return results, fig

    def save_plot(self, filepath: str, dpi: int = 150, **kwargs):
        if self.figure is None:
            self.plot()
        self.figure.savefig(filepath, dpi=dpi, bbox_inches='tight', **kwargs)

    def get_summary(self) -> str:
        if self.results is None:
            return "No analysis results available. Run analyze() first."
        return self.analyzer.get_results_summary()

    def close_plot(self):
        if self.figure is not None:
            plt.close(self.figure)
            self.figure = None
