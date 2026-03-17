"""Image Processing utilities used across Alexandria analyzers.

This module collects lightweight, well-documented image processing
helpers commonly needed by CatPhan analyzers and plotters.
"""

import numpy as np
from scipy import ndimage
from scipy.interpolate import interpn
from typing import Tuple, Optional


class ImageProcessor:
    @staticmethod
    def apply_gaussian_filter(image: np.ndarray, sigma: float = 1.0) -> np.ndarray:
        return ndimage.gaussian_filter(image, sigma=sigma)

    @staticmethod
    def extract_profile(image: np.ndarray, start: Tuple[float, float], end: Tuple[float, float], n_points: int = 100) -> np.ndarray:
        x_samples = np.linspace(start[0], end[0], n_points)
        y_samples = np.linspace(start[1], end[1], n_points)
        ny, nx = image.shape
        x_grid = np.linspace(0, nx - 1, nx)
        y_grid = np.linspace(0, ny - 1, ny)
        sample_coords = np.vstack((y_samples, x_samples)).T
        vals = interpn((y_grid, x_grid), image, sample_coords, method='linear', bounds_error=False, fill_value=0.0)
        return vals

    @staticmethod
    def threshold_image(image: np.ndarray, threshold: float, mode: str = 'above') -> np.ndarray:
        if mode == 'above':
            return image > threshold
        elif mode == 'below':
            return image < threshold
        else:
            raise ValueError("mode must be 'above' or 'below'")

    @staticmethod
    def estimate_noise(image: np.ndarray, roi_center: Optional[Tuple[float, float]] = None, roi_size: int = 50) -> float:
        ny, nx = image.shape
        if roi_center is None:
            cx, cy = nx // 2, ny // 2
        else:
            cx, cy = roi_center
        half_size = roi_size // 2
        x_start = max(0, int(cx - half_size))
        x_end = min(nx, int(cx + half_size))
        y_start = max(0, int(cy - half_size))
        y_end = min(ny, int(cy + half_size))
        roi = image[y_start:y_end, x_start:x_end]
        return float(np.std(roi))

    @staticmethod
    def find_edges(image: np.ndarray, method: str = 'sobel') -> np.ndarray:
        if method == 'sobel':
            sx = ndimage.sobel(image, axis=0)
            sy = ndimage.sobel(image, axis=1)
        elif method == 'prewitt':
            sx = ndimage.prewitt(image, axis=0)
            sy = ndimage.prewitt(image, axis=1)
        elif method == 'scharr':
            kx = np.array([[-3, 0, 3], [-10, 0, 10], [-3, 0, 3]])
            ky = kx.T
            sx = ndimage.convolve(image, kx)
            sy = ndimage.convolve(image, ky)
        else:
            raise ValueError("method must be 'sobel', 'prewitt', or 'scharr'")
        return np.hypot(sx, sy)
