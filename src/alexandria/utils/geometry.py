"""
Geometry Utilities for CatPhan Phantom Analysis

Provides functions for finding phantom centers, rotations, and geometric measurements.
"""

import numpy as np
from scipy.interpolate import interpn
from scipy.signal import find_peaks
from typing import Tuple, List, Optional, Union


class CatPhanGeometry:
    @staticmethod
    def find_center(image: np.ndarray, threshold: float = 400) -> Tuple[List[float], List[np.ndarray]]:
        sz       = np.array(image.shape)
        matrix_c = (np.round(sz[0]/2), np.round(sz[1]/2))
        px = image[int(matrix_c[0]), :]
        py = image[:, int(matrix_c[1])]
        offset = 1
        try:
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
        except StopIteration:
            threshold = 300
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
        szx = x2 - x1
        szy = y2 - y1
        center = [(x1 + x2) / 2, (y1 + y2) / 2]
        outer_r = (szx + szy) / 4
        t = np.linspace(0, 2*np.pi, 100)
        outer_x = outer_r * np.cos(t) + center[0]
        outer_y = outer_r * np.sin(t) + center[1]
        return center, [outer_x, outer_y]

    @staticmethod
    def select_optimal_ctp528_slices(dicom_set: List, target_index: int, search_range: int = 2) -> Tuple[np.ndarray, np.ndarray, float]:
        return

    @staticmethod
    def calculate_slice_thickness(image: np.ndarray, pixel_spacing: float, center: Tuple[float, float]) -> float:
        return 5.0


def circular_roi_mask(shape: Tuple[int, int], center: Tuple[float, float], radius: float) -> np.ndarray:
    ny, nx = shape
    Y, X = np.ogrid[:ny, :nx]
    cx, cy = center
    dist = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    return dist <= radius


def compute_phantom_boundary(image: np.ndarray, center: Tuple[float, float], pixel_spacing: Optional[float] = None, threshold: float = -900, fallback_threshold: float = -900) -> Tuple[Tuple[float, float], Tuple[np.ndarray, np.ndarray]]:
    if image is None or center is None:
        return ((0, 0), (np.array([]), np.array([])))
    center_row = int(round(center[1]))
    center_col = int(round(center[0]))
    center_row = max(0, min(center_row, image.shape[0] - 1))
    center_col = max(0, min(center_col, image.shape[1] - 1))
    px = image[center_row, :]
    py = image[:, center_col]
    offset = 1
    try:
        x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
        y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
        x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
        y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
    except StopIteration:
        try:
            threshold = fallback_threshold
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
        except StopIteration:
            if pixel_spacing:
                radius_px = 100.0 / pixel_spacing
                t = np.linspace(0, 2*np.pi, 100)
                outer_x = radius_px * np.cos(t) + center[0]
                outer_y = radius_px * np.sin(t) + center[1]
                return (center, (outer_x, outer_y))
            return ((0, 0), (np.array([]), np.array([])))
    szx = x2 - x1
    szy = y2 - y1
    detected_center = ((x1 + x2) / 2, (y1 + y2) / 2)
    outer_r = (szx + szy) / 4
    t = np.linspace(0, 2*np.pi, 100)
    outer_x = outer_r * np.cos(t) + detected_center[0]
    outer_y = outer_r * np.sin(t) + detected_center[1]
    return (detected_center, (outer_x, outer_y))


def draw_boundary(center: Tuple[float, float], diameter_x_px: Optional[float], diameter_y_px: Optional[float], n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    if diameter_x_px is None or diameter_y_px is None:
        return np.array([]), np.array([])
    rx = diameter_x_px / 2
    ry = diameter_y_px / 2
    t = np.linspace(0, 2 * np.pi, n_points)
    x_coords = rx * np.cos(t) + center[0]
    y_coords = ry * np.sin(t) + center[1]
    return x_coords, y_coords


def find_center_edge_detection(img: np.ndarray, threshold: float = -900, fallback_threshold: float = -900, return_diameters: bool = False):
    sz = np.array(img.shape)
    matrix_c = (int(np.round(sz[0]/2)), int(np.round(sz[1]/2)))
    px = img[matrix_c[0], :]
    py = img[:, matrix_c[1]]
    offset = 1
    try:
        x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
        y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
        x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
        y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
    except StopIteration:
        threshold = fallback_threshold
        try:
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
        except StopIteration:
            if return_diameters:
                return matrix_c[0], matrix_c[1], None, None
            return matrix_c[0], matrix_c[1]
    center_col = int((x1 + x2) / 2)
    center_row = int((y1 + y2) / 2)
    if return_diameters:
        diameter_x = float(x2 - x1)
        diameter_y = float(y2 - y1)
        return center_row, center_col, diameter_y, diameter_x
    return center_row, center_col


def find_rotation(image: np.ndarray, center: Optional[Tuple[float, float]], pixel_spacing: Union[float, Tuple[float, float], List[float]], insert_radius_mm: float = 58.5, edge_threshold: float = 100.0, center_threshold: float = 30, iterations: int = 5, profile_length: int = 25, granularity: int = 3, interp_kwargs: Optional[dict] = None, initial_angle_deg: float = 0.0):
    if interp_kwargs is None:
        interp_kwargs = {'bounds_error': False, 'fill_value': 0}
    if isinstance(pixel_spacing, (list, tuple, np.ndarray)):
        space = float(pixel_spacing[0])
    else:
        space = float(pixel_spacing)
    if center is None:
        try:
            r_row, r_col = find_center_edge_detection(image)
            center = (float(r_col), float(r_row))
        except Exception:
            raise ValueError("Center must be provided or detectable via edge detection")
    h_img, w_img = image.shape[:2]
    ring_r = insert_radius_mm / space
    _p90 = (ring_r * np.cos(np.radians(90 + initial_angle_deg)) + center[0], ring_r * np.sin(np.radians(90 + initial_angle_deg)) + center[1])
    _p270 = (ring_r * np.cos(np.radians(270 + initial_angle_deg)) + center[0], ring_r * np.sin(np.radians(270 + initial_angle_deg)) + center[1])
    ct = _p270
    cb = _p90
    x = np.linspace(0, h_img - 1, h_img)
    y = np.linspace(0, w_img - 1, w_img)
    def _find_insert_center(roi_pos):
        x_horiz = np.linspace(roi_pos[0] - profile_length, roi_pos[0] + profile_length, profile_length * granularity)
        x_vert = np.linspace(roi_pos[1] - profile_length, roi_pos[1] + profile_length, profile_length * granularity)
        prof_h = np.zeros(len(x_horiz))
        prof_v = np.zeros(len(x_vert))
        for i in range(len(x_horiz)):
            prof_h[i] = interpn((x, y), image, [roi_pos[1], x_horiz[i]], **interp_kwargs)
        for i in range(len(x_vert)):
            prof_v[i] = interpn((x, y), image, [x_vert[i], roi_pos[0]], **interp_kwargs)
        dh = np.diff(prof_h)
        dv = np.diff(prof_v)
        peaks_h, _ = find_peaks(np.abs(dh), height=edge_threshold)
        peaks_v, _ = find_peaks(np.abs(dv), height=edge_threshold)
        if len(peaks_h) >= 2 and len(peaks_v) >= 2:
            offset_len = len(x_horiz) / 2
            mid_h = np.mean(peaks_h) - offset_len
            mid_v = np.mean(peaks_v) - offset_len
            return (roi_pos[0] + mid_h, roi_pos[1] + mid_v)
        return roi_pos
    ct_orig = ct
    cb_orig = cb
    ct_old = ct
    cb_old = cb
    for _ in range(iterations):
        try:
            ct_new = _find_insert_center(ct)
            cb_new = _find_insert_center(cb)
        except Exception:
            ct, cb = ct_orig, cb_orig
            break
        if (abs(ct_new[0] - ct_old[0]) > center_threshold or abs(ct_new[1] - ct_old[1]) > center_threshold or abs(cb_new[0] - cb_old[0]) > center_threshold or abs(cb_new[1] - cb_old[1]) > center_threshold):
            ct, cb = ct_orig, cb_orig
            break
        ct_old, cb_old = ct_new, cb_new
        ct, cb = ct_new, cb_new
    tx = ct[0] - cb[0]
    ty = ct[1] - cb[1]
    rotation_angle = np.degrees(np.arctan2(-ty, tx))
    rotation_from_y = float(rotation_angle) - 90.0
    return rotation_from_y, ct, cb
