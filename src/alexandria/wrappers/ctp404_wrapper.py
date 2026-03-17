"""
CTP404 Module Wrapper

Combines CTP404Analyzer with CTP404Plotter for convenient
9-ROI sensitometry analysis and visualization in a single class.
"""

import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict, Any
from ..analyzers.ctp404 import CTP404Analyzer


class CTP404ModuleReporter:
    def __init__(
        self,
        image: Optional[Any] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,
        dicom_set: Optional[List[Any]] = None,
        slice_index: Optional[int] = None,
        rotation_offset: float = 0.0,
        roi_radius: float = 3.5,
        material_distance: float = 58.5,
    ):
        self.analyzer = CTP404Analyzer(
            image=image,
            center=center,
            pixel_spacing=pixel_spacing,
            spacing=spacing,
            dicom_set=dicom_set,
            slice_index=slice_index,
            rotation_offset=rotation_offset,
            roi_radius=roi_radius,
            material_distance=material_distance,
        )
        self.results = None
        self.figure = None

    def analyze(self, verbose: bool = True) -> Dict[str, Any]:
        self.results = self.analyzer.analyze(verbose=verbose)
        return self.results

    def plot(self, show: bool = False, **kwargs) -> plt.Figure:
        if self.results is None:
            self.analyze(verbose=False)
        from ..plotters.ctp404_plotter import CTP404Plotter
        plotter = CTP404Plotter(self.analyzer)
        self.figure = plotter.plot(**kwargs)
        if show:
            plt.show()
        return self.figure

    def analyze_and_plot(self, verbose: bool = True, show: bool = False, **kwargs) -> Tuple[Dict[str, Any], plt.Figure]:
        results = self.analyze(verbose=verbose)
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
