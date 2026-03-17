"""
Detailed Uniformity Analyzer

Samples pixel values along concentric circular profiles and records
angle/value pairs for each radius.
"""

from typing import Any, Dict, List, Optional, Tuple, Callable

import numpy as np
from scipy.ndimage import map_coordinates


class DetailedUniformityAnalyzer:
    """
    Analyzer that samples concentric circular profiles for uniformity checks.
    """
    def __init__(
        self,
        image: Optional[np.ndarray] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,
        dicom_set: Optional[List[Any]] = None,
        slice_index: Optional[int] = None,
        radii_mm: Optional[List[float]] = None,
        sample_step_mm: float = 1.0,
        n_samples: int = 360,
        center_finder: Optional[Callable[..., Tuple]] = None,
        center_finder_kwargs: Optional[Dict[str, Any]] = None,
        center_threshold: float = 400.0,
        center_threshold_fallback: float = -900.0,
    ):
        self.dicom_mode = dicom_set is not None and slice_index is not None
        self.dicom_set = dicom_set
        self.slice_index = slice_index
        self.averaged_image = None

        self.image = np.array(image, dtype=float) if image is not None else None
        self.center = center

        if pixel_spacing is not None:
            self.pixel_spacing = float(pixel_spacing)
        elif spacing is not None:
            self.pixel_spacing = float(spacing)
        else:
            self.pixel_spacing = None

        self.radii_mm = radii_mm if radii_mm is not None else [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0, 60.0, 65.0, 70.0]
        self.sample_step_mm = float(sample_step_mm)
        self.n_samples = int(n_samples)

        self.center_finder = center_finder
        self.center_finder_kwargs = center_finder_kwargs or {}
        self.center_threshold = center_threshold
        self.center_threshold_fallback = center_threshold_fallback

        self.profile_data: List[Dict[str, Any]] = []
        self.results: Dict[str, Any] = {}

    def prepare_image(self) -> np.ndarray:
        """
        Prepare image for analysis.
        """
        if self.dicom_mode:
            idx = self.slice_index
            if idx is None:
                raise ValueError("slice_index must be provided for DICOM mode")

            if idx <= 0:
                im1 = self.dicom_set[idx].pixel_array
                im2 = self.dicom_set[idx + 1].pixel_array
                im3 = self.dicom_set[idx + 1].pixel_array
            elif idx >= len(self.dicom_set) - 1:
                im1 = self.dicom_set[idx].pixel_array
                im2 = self.dicom_set[idx - 1].pixel_array
                im3 = self.dicom_set[idx - 1].pixel_array
            else:
                im1 = self.dicom_set[idx].pixel_array
                im2 = self.dicom_set[idx + 1].pixel_array
                im3 = self.dicom_set[idx - 1].pixel_array

            self.averaged_image = (im1 + im2 + im3) / 3
            self.image = self.averaged_image

            if self.pixel_spacing is None:
                self.pixel_spacing = float(self.dicom_set[idx].PixelSpacing[0])

            return self.averaged_image

        if self.image is None:
            raise ValueError("Image must be provided in single-image mode")

        return self.image

    def _compute_center(self) -> Tuple[float, float, Optional[float], Optional[float]]:
        if self.center is not None:
            return float(self.center[0]), float(self.center[1]), None, None

        from alexandria.utils import find_center_edge_detection

        def _unpack(value: Any) -> Tuple[int, int, Optional[float], Optional[float]]:
            if isinstance(value, (tuple, list)):
                if len(value) >= 4:
                    return int(value[0]), int(value[1]), value[2], value[3]
                if len(value) == 2:
                    return int(value[0]), int(value[1]), None, None
            raise ValueError(
                "center_finder must return (row, col) or (row, col, diameter_y_px, diameter_x_px)"
            )

        if self.center_finder is not None:
            result = self.center_finder(self.image, **self.center_finder_kwargs)
            row, col, diameter_y, diameter_x = _unpack(result)
        else:
            row, col, diameter_y, diameter_x = find_center_edge_detection(
                self.image,
                threshold=self.center_threshold,
                fallback_threshold=self.center_threshold_fallback,
                return_diameters=True
            )

        self.center = (float(col), float(row))
        return float(col), float(row), diameter_y, diameter_x

    def _resolve_radii_mm(self) -> List[float]:
        max_radius = float(max(self.radii_mm)) if self.radii_mm else 0.0
        if max_radius <= 0.0:
            return []
        count = int(round(max_radius / self.sample_step_mm))
        return [self.sample_step_mm * idx for idx in range(1, count + 1)]

    def _plot_radii_mm(self) -> List[float]:
        return [float(r) for r in self.radii_mm]

    def _sample_circle(self, center_x: float, center_y: float, radius_px: float) -> Tuple[np.ndarray, np.ndarray]:
        angles_deg = np.linspace(0.0, 360.0, self.n_samples, endpoint=False)
        angles_rad = np.deg2rad(angles_deg)

        xs = center_x + radius_px * np.cos(angles_rad)
        ys = center_y + radius_px * np.sin(angles_rad)

        coords = np.vstack([ys, xs])
        values = map_coordinates(self.image, coords, order=1, mode='nearest')
        return angles_deg, values

    def analyze(self) -> Dict[str, Any]:
        if self.image is None:
            self.prepare_image()

        if self.pixel_spacing is None:
            raise ValueError("Pixel spacing must be provided")

        center_x, center_y, _, _ = self._compute_center()
        radii_mm = self._resolve_radii_mm()

        self.profile_data = []
        profile_results = []

        for radius_mm in radii_mm:
            radius_px = radius_mm / self.pixel_spacing
            angles_deg, values = self._sample_circle(center_x, center_y, radius_px)

            profile = {
                "radius_mm": float(radius_mm),
                "radius_px": float(radius_px),
                "angles_deg": angles_deg,
                "values": values,
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }
            self.profile_data.append(profile)
            profile_results.append({
                "radius_mm": profile["radius_mm"],
                "radius_px": profile["radius_px"],
                "angles_deg": profile["angles_deg"].tolist(),
                "values": profile["values"].tolist(),
                "mean": profile["mean"],
                "std": profile["std"],
                "min": profile["min"],
                "max": profile["max"],
            })

        self.results = {
            "center": [center_x, center_y],
            "pixel_spacing": float(self.pixel_spacing),
            "n_samples": self.n_samples,
            "profiles": profile_results,
            "plot_radii_mm": self._plot_radii_mm(),
            "sample_step_mm": float(self.sample_step_mm),
        }

        return self.results
