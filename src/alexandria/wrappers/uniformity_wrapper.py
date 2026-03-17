"""
Uniformity Module Wrapper

Combines UniformityAnalyzer with UniformityPlotter for convenient
analysis and visualization in a single class.
"""

import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict, Any
from ..analyzers.uniformity import UniformityAnalyzer


class UniformityModuleReporter:
    def __init__(
        self,
        image: Optional[Any] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,
        dicom_set: Optional[List] = None,
        slice_index: Optional[int] = None,
        roi_box_size: float = 15.0,
        roi_offset: float = 50.0,
    ):
        self.analyzer = UniformityAnalyzer(
            image=image,
            center=center,
            pixel_spacing=pixel_spacing,
            spacing=spacing,
            dicom_set=dicom_set,
            slice_index=slice_index,
            roi_box_size=roi_box_size,
            roi_offset=roi_offset,
        )
        self.results = None
        self.figure = None

    def analyze(self, verbose: bool = True) -> Dict[str, Any]:
        self.results = self.analyzer.analyze(verbose=verbose)
        return self.results

    def plot(self, show: bool = False, **kwargs) -> plt.Figure:
        if self.results is None:
            self.analyze(verbose=False)
        from ..plotters.uniformity_plotter import UniformityPlotter
        plotter = UniformityPlotter(self.analyzer)
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
