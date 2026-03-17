"""
Unified High Contrast Analyzer for CatPhan Phantom Analysis

This module combines the high-contrast/resolution analysis functionality from both
catphan404 and XVI-CatPhan projects, providing a comprehensive analyzer for CTP528
line pair module with both single-image and DICOM-series modes.
"""

import numpy as np
from scipy.interpolate import interpn, interp1d
from scipy.signal import find_peaks
from typing import Optional, Tuple, List, Dict, Any


class HighContrastAnalyzer:
    """
    Analyzer for CatPhan's CTP528 high-contrast (line pair) module.
    """

    def __init__(
        self,
        image: Optional[np.ndarray] = None,
        pixel_spacing: Optional[float] = None,
        center: Optional[Tuple[float, float]] = None,
        t_offset_deg: float = 0.0,
        rotation_offset: Optional[float] = None,  # Alias for t_offset_deg
        dicom_set: Optional[List] = None,
        slice_index: Optional[int] = None,
        lp_r_mm: float = 48.0,
        samples_per_segment: int = 50,
        center_threshold: float = -980,
        center_threshold_fallback: float = -900.0,
    ):
        # Mode detection
        self.dicom_mode = dicom_set is not None and slice_index is not None
        self.dicom_set = dicom_set
        self.slice_index = slice_index
        self.averaged_image = None

        if image is not None:
            self.image = np.array(image, dtype=float)
        else:
            self.image = None

        self.pixel_spacing             = float(pixel_spacing) if pixel_spacing is not None else None
        self.center                    = center
        self.center_threshold          = center_threshold
        self.center_threshold_fallback = center_threshold_fallback

        if rotation_offset is not None:
            self.t_offset_deg = float(rotation_offset)
        else:
            self.t_offset_deg = float(t_offset_deg)

        self.lp_r_mm = float(lp_r_mm)
        self.samples_per_segment = int(samples_per_segment)

        if center is not None:
            self.center_x = float(center[0])
            self.center_y = float(center[1])
        else:
            self.center_x = None
            self.center_y = None

        self.theta_deg = np.array([10, 40, 62, 85, 103, 121, 140, 157, 173, 186]) + self.t_offset_deg
        self.npeaks = [[1, 2], [2, 3], [3, 4], [4, 4], [5, 4], [6, 5], [7, 5], [8, 5], [9, 5], [10, 5]]

        self.lpx = None
        self.lpy = None

        self.per_pair_mtf = []
        self.profiles = []
        self.peaks_max = []
        self.peaks_min = []
        self.peaks_combined = []
        self.lp_x = []
        self.lp_y = []
        self.nMTF = None
        self.lp_axis = None
        self.mtf_points = {}
        self.line_pair_profiles = None

        self.results = {}
        self.mtf = None
        self.lp_frequencies = None

    def prepare_image(self):
        if self.dicom_mode:
            try:
                from ..utils.geometry import CatPhanGeometry

                im, means, z_mean = CatPhanGeometry.select_optimal_ctp528_slices(
                    self.dicom_set, self.slice_index
                )
                self.averaged_image = im
                self.image = im
            except Exception:
                idx = self.slice_index
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

                self.averaged_image = (im1.astype(float) + im2.astype(float) + im3.astype(float)) / 3
                self.image = self.averaged_image

            if self.pixel_spacing is None:
                self.pixel_spacing = float(self.dicom_set[self.slice_index].PixelSpacing[0])

            if self.center is not None:
                self.center_x = float(self.center[0])
                self.center_y = float(self.center[1])

            return self.averaged_image
        else:
            if self.image is None:
                raise ValueError("Image must be provided in single-image mode")
            return self.image

    def _compute_centers(self):
        if self.center_x is None or self.center_y is None:
            raise ValueError("Center coordinates must be provided")
        if self.pixel_spacing is None:
            raise ValueError("Pixel spacing must be provided")

        r_pixels = self.lp_r_mm / self.pixel_spacing
        thetas_rad = np.deg2rad(self.theta_deg)
        self.lpx = r_pixels * np.cos(thetas_rad) + self.center_x
        self.lpy = r_pixels * np.sin(thetas_rad) + self.center_y

    def _get_MTF_for_pair(self, x_coords: Tuple[float, float], y_coords: Tuple[float, float], npeaks_expected: List[int]):
        x1 = np.linspace(x_coords[0], x_coords[1], self.samples_per_segment)
        y1 = np.linspace(y_coords[0], y_coords[1], self.samples_per_segment)
        f1 = np.zeros(len(x1))

        ny, nx = self.image.shape
        x = np.linspace(0, nx - 1, nx)
        y = np.linspace(0, ny - 1, ny)

        for i in range(len(x1)):
            f1[i] = interpn((y, x), self.image, [[y1[i], x1[i]]], method='linear', bounds_error=False, fill_value=0.0)[0]

        df1 = np.diff(f1)

        h = 50
        peaks_max1, _ = find_peaks(df1, height=h)
        peaks_min1, _ = find_peaks(-df1, height=h)

        while (len(peaks_max1) < npeaks_expected[1]) or (len(peaks_min1) < npeaks_expected[1]):
            if h <= 10:
                return 0.0, f1, peaks_max1, peaks_min1, np.array([], dtype=int), x_coords, y_coords
            h -= 1
            peaks_max1, _ = find_peaks(df1, height=h)
            peaks_min1, _ = find_peaks(-df1, height=h)

        peaks1 = np.hstack((peaks_max1, peaks_min1))
        peaks1 = np.array(sorted(peaks1))

        idxmax = []
        idxmin = []
        Imax = []
        Imin = []
        offset = 1

        for k in range(len(peaks1) - 1):
            if k % 2 == 0:
                tmp_idx = np.array(f1[peaks1[k] - offset:peaks1[k + 1] + offset]).argmax()
                idx_at = tmp_idx - offset + peaks1[k]
                idxmax.append(idx_at)
                Imax.append(f1[idx_at])
            else:
                tmp_idx = np.array(f1[peaks1[k] - offset:peaks1[k + 1] + offset]).argmin()
                idx_at = tmp_idx - offset + peaks1[k]
                idxmin.append(idx_at)
                Imin.append(f1[idx_at])

        if (len(Imax) == 0) or (len(Imin) == 0):
            MTF_value = 0.0
        else:
            MTF_value = (np.mean(Imax) - np.mean(Imin)) / (np.mean(Imax) + np.mean(Imin))

        return float(MTF_value), f1, peaks_max1, peaks_min1, peaks1, x_coords, y_coords

    def analyze(self, write_log: bool = False, verbose: bool = True) -> Dict[str, Any]:
        if self.image is None:
            self.prepare_image()

        if self.center is None:
            from alexandria.utils import find_center_edge_detection, compute_phantom_boundary, draw_boundary
            center_row, center_col, diameter_y, diameter_x = find_center_edge_detection(
                self.image,
                threshold=self.center_threshold,
                fallback_threshold=self.center_threshold_fallback,
                return_diameters=True
            )
            self.center = (center_col, center_row)
            self.center_x = float(self.center[0])
            self.center_y = float(self.center[1])

            boundary_x, boundary_y = draw_boundary(self.center, diameter_x, diameter_y)
            if len(boundary_x) == 0:
                _, (boundary_x, boundary_y) = compute_phantom_boundary(
                    self.image,
                    self.center,
                    self.pixel_spacing,
                    threshold=self.center_threshold,
                    fallback_threshold=self.center_threshold_fallback
                )
            self.boundary = {
                'x': boundary_x.tolist() if len(boundary_x) > 0 else [],
                'y': boundary_y.tolist() if len(boundary_y) > 0 else []
            }

        self._compute_centers()

        n_pairs = len(self.lpx) - 1
        per_pair_mtf = []
        profiles = []
        pmax_list = []
        pmin_list = []
        pcomb_list = []
        lp_x_list = []
        lp_y_list = []

        for i in range(n_pairs):
            npeaks = self.npeaks[i]
            x_coords = (self.lpx[i], self.lpx[i + 1])
            y_coords = (self.lpy[i], self.lpy[i + 1])

            mtf_val, f1, pmax, pmin, pcomb, tmpx, tmpy = self._get_MTF_for_pair(
                x_coords, y_coords, npeaks
            )
            per_pair_mtf.append(mtf_val)
            profiles.append(f1)
            pmax_list.append(pmax)
            pmin_list.append(pmin)
            pcomb_list.append(pcomb)
            lp_x_list.append(tmpx)
            lp_y_list.append(tmpy)

        self.per_pair_mtf = per_pair_mtf
        self.profiles = profiles
        self.line_pair_profiles = profiles
        self.peaks_max = pmax_list
        self.peaks_min = pmin_list
        self.peaks_combined = pcomb_list
        self.lp_x = lp_x_list
        self.lp_y = lp_y_list

        mtf_array = np.array(self.per_pair_mtf, dtype=float)
        max_val = mtf_array.max() if mtf_array.size > 0 else 0.0
        if max_val > 0:
            self.nMTF = mtf_array / max_val
        else:
            self.nMTF = mtf_array.copy()

        if len(self.per_pair_mtf) > 0:
            self.lp_axis = np.linspace(1, len(self.per_pair_mtf), len(self.per_pair_mtf)) / 10.0
        else:
            self.lp_axis = np.array([])

        if self.nMTF.size > 0 and self.lp_axis.size > 0:
            targets = (0.8, 0.5, 0.3, 0.1)
            try:
                fMTF = np.interp(targets, self.nMTF[::-1], self.lp_axis[::-1])
                self.mtf_points = {
                    "MTF80": float(fMTF[0]),
                    "MTF50": float(fMTF[1]),
                    "MTF30": float(fMTF[2]),
                    "MTF10": float(fMTF[3]),
                }
            except Exception:
                self.mtf_points = {"MTF80": np.nan, "MTF50": np.nan, "MTF30": np.nan, "MTF10": np.nan}
        else:
            self.mtf_points = {"MTF80": np.nan, "MTF50": np.nan, "MTF30": np.nan, "MTF10": np.nan}

        self.results = {
            'mtf_80': self.mtf_points["MTF80"],
            'mtf_50': self.mtf_points["MTF50"],
            'mtf_30': self.mtf_points["MTF30"],
            'mtf_10': self.mtf_points["MTF10"],
            'mtf_array': self.nMTF,
            'lp_frequencies': self.lp_axis,
            'line_pair_x': self.lpx,
            'line_pair_y': self.lpy
        }
        
        self.mtf = self.nMTF
        self.lp_frequencies = self.lp_axis

        return self.results

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pixel_spacing_mm": self.pixel_spacing,
            "image_shape": self.image.shape if self.image is not None else None,
            "centers_x": self.lpx.tolist() if self.lpx is not None else None,
            "centers_y": self.lpy.tolist() if self.lpy is not None else None,
            "per_pair_mtf": [float(x) for x in self.per_pair_mtf],
            "MTF10_lp_per_mm": self.mtf_points.get("MTF10"),
            "MTF30_lp_per_mm": self.mtf_points.get("MTF30"),
            "MTF50_lp_per_mm": self.mtf_points.get("MTF50"),
            "MTF80_lp_per_mm": self.mtf_points.get("MTF80"),
            "lp_mm": self.lp_axis.tolist() if self.lp_axis is not None else None,
            "nmtf": self.nMTF.tolist() if self.nMTF is not None else None,
        }
