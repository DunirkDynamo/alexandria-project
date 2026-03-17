"""
Unified CTP515 Low-Contrast Analyzer for CatPhan Phantom Analysis

This module provides low-contrast detectability analysis for the CatPhan
CTP515 module. Analyzes circular inserts of varying diameters (15, 9, 8, 7, 6, 5 mm)
to measure Contrast-to-Noise Ratio (CNR) and contrast detectability.
"""

# Sphinx/autodoc: expanded docstrings and inline comments for clarity.

import numpy as np
import math
from typing import Optional, Tuple, List, Dict, Any, Callable


class CTP515Analyzer:
    """
    Low-contrast detectability analyzer for CatPhan CTP515 module.

    Detects and analyzes six low-contrast circular inserts of varying diameters
    positioned at fixed angles and distance from center. Computes Contrast-to-Noise
    Ratio (CNR) and contrast percentage for each ROI relative to a background region.

    CNR quantifies detectability: higher values indicate the insert is more
    easily distinguished from background noise.

    Supports both single-image mode and DICOM-series mode with 3-slice averaging.

    Attributes:
        image (np.ndarray): 2D CT image of the low-contrast module.
        center (tuple): (x, y) center of phantom in pixels.
        pixel_spacing (float): Pixel spacing in mm.
        angle_offset (float): Angular offset for ROI positioning in degrees.
        results (dict): Analysis results populated by analyze().
    """

    # ROI specifications (angles in degrees, sizes in mm)
    ROI_ANGLES = [
        -87.4,
        -69.1,#-105.7,
        -52.7,
        -38.5,
        -25.1,
        -12.9,
    ]

    ROI_DISTANCE_MM = 50  # Distance from center to ROI centers
    ROI_RADII_MM = [6, 3.5, 3, 2.5, 2, 1.5]  # ROI radii for 15, 9, 8, 7, 6, 5 mm diameters

    # ROI settings mapped by diameter
    ROI_SETTINGS = {
        "15": {"angle_idx": 0, "radius_mm": 6},
        "9": {"angle_idx": 1, "radius_mm": 3.5},
        "8": {"angle_idx": 2, "radius_mm": 3},
        "7": {"angle_idx": 3, "radius_mm": 2.5},
        "6": {"angle_idx": 4, "radius_mm": 2},
        "5": {"angle_idx": 5, "radius_mm": 1.5},
    }

    def __init__(
        self,
        image: Optional[np.ndarray] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        angle_offset: float = 0.0,
        dicom_set: Optional[List] = None,
        slice_index: Optional[int] = None,
        center_finder: Optional[Callable[..., Tuple]] = None,
        center_finder_kwargs: Optional[Dict[str, Any]] = None,
        center_threshold: float = 400.0,
        center_threshold_fallback: float = -900.0,
    ):
        """
        Initialize the CTP515 analyzer.

        Supports two initialization modes:
        1. Single-image mode: Provide image, center, and pixel_spacing
        2. DICOM-series mode: Provide dicom_set and slice_index

        Args:
            image: 2D CT image of the low-contrast module (single-image mode).
            center: (x, y) coordinates of phantom center in pixels.
            pixel_spacing: Pixel spacing in mm.
            angle_offset: Angular offset for ROI positioning in degrees.
            dicom_set: List of DICOM dataset objects (DICOM-series mode).
            slice_index: Index of the CTP515 slice in the dataset.
            center_finder: Optional callable to compute center from image.
                           Should return (row, col) or (row, col, diameter_y_px, diameter_x_px).
            center_finder_kwargs: Optional dict of kwargs to pass to center_finder.
            center_threshold: Threshold used by the default center finder.
            center_threshold_fallback: Fallback threshold if the primary fails.
        """
        # Mode detection
        self.dicom_mode = dicom_set is not None and slice_index is not None

        # DICOM-series mode attributes
        self.dicom_set = dicom_set
        self.slice_index = slice_index
        self.averaged_image = None

        # Image and parameters
        if image is not None:
            self.image = np.array(image, dtype=float)
        else:
            self.image = None

        self.center = center
        self.pixel_spacing = pixel_spacing
        self.angle_offset = angle_offset
        self.center_finder = center_finder
        self.center_finder_kwargs = center_finder_kwargs or {}
        self.center_threshold = center_threshold
        self.center_threshold_fallback = center_threshold_fallback

        # Results storage
        self.results = {}

    def prepare_image(self):
        """
        Prepare image for analysis.

        In DICOM mode: Create 3-slice averaged image for improved SNR.
        In single-image mode: Use the provided image directly.

        Returns:
            Prepared image array
        """
        if self.dicom_mode:
            # DICOM-series mode: 3-slice averaging. Use neighbor repetition at
            # series boundaries to avoid IndexError; convert to float for safety.
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

            self.averaged_image = (im1.astype(float) + im2.astype(float) + im3.astype(float)) / 3.0
            self.image = self.averaged_image

            # Extract pixel spacing from DICOM if not provided
            if self.pixel_spacing is None:
                self.pixel_spacing = float(self.dicom_set[self.slice_index].PixelSpacing[0])

            return self.averaged_image
        else:
            # Single-image mode
            if self.image is None:
                raise ValueError("Image must be provided in single-image mode")
            return self.image

    def _validate_inputs(self):
        """
        Validate the input parameters to prevent runtime errors.

        Raises:
            TypeError: If inputs are not of the correct type.
            ValueError: If image is not 2D.
        """
        if not isinstance(self.image, np.ndarray):
            raise TypeError("image must be a numpy.ndarray")
        if self.image.ndim != 2:
            raise ValueError("image must be a 2D array")
        if not (isinstance(self.center, tuple) and len(self.center) == 2):
            raise TypeError("center must be a tuple of (x, y)")
        if not isinstance(self.pixel_spacing, (float, int)):
            raise TypeError("pixel_spacing must be a float or int")
        if not isinstance(self.angle_offset, (float, int)):
            raise TypeError("angle_offset must be a float or int")

    def _circular_roi_mask(
        self,
        shape: Tuple[int, int],
        center: Tuple[float, float],
        radius: float
    ) -> np.ndarray:
        """
        Create a circular boolean mask.

        Args:
            shape: Image shape (height, width).
            center: (x, y) coordinates of circle center.
            radius: Circle radius in pixels.

        Returns:
            Boolean mask array.
        """
        ny, nx = shape
        y, x = np.ogrid[:ny, :nx]
        cx, cy = center
        dist_from_center = np.sqrt((x - cx)**2 + (y - cy)**2)
        return dist_from_center <= radius

    def analyze(self, verbose: bool = True) -> Dict[str, Any]:
        """
        Perform low-contrast detectability analysis.

        This method:
        1. Defines ROI locations based on predefined angles and distances.
        2. For each ROI, creates a circular mask and computes mean/std.
        3. Computes a common background ROI for noise reference.
        4. Calculates CNR for each ROI against the background.
        5. Returns a summary of detected ROIs and their metrics.

        Args:
            verbose: Whether to print progress information.

        Returns:
            Dict: Contains 'n_detected' (int) and 'blobs' (dict of blob stats).
                  Each blob entry has position, size, means, std, and CNR.
        """
        # Prepare image if needed
        if self.image is None:
            self.prepare_image()

        # Compute center if not provided
        if self.center is None:
            from alexandria.utils import find_center_edge_detection, compute_phantom_boundary, draw_boundary

            def _unpack_center_result(value: Any) -> Tuple[int, int, Optional[float], Optional[float]]:
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
                center_row, center_col, diameter_y, diameter_x = _unpack_center_result(result)
            else:
                center_row, center_col, diameter_y, diameter_x = find_center_edge_detection(
                    self.image,
                    threshold=self.center_threshold,
                    fallback_threshold=self.center_threshold_fallback,
                    return_diameters=True
                )
            self.center = (center_col, center_row)  # Convert to (x, y)

            # Draw boundary from edge-derived diameters when available
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

        # Validate inputs
        self._validate_inputs()

        # Get image dimensions and center coordinates
        # center is passed as (x, y) = (col, row)
        ny, nx = self.image.shape
        cx, cy = self.center[0], self.center[1]

        # Initialize results dictionary
        results = {}

        # Compute background statistics from a common background ROI
        # Background ROI: 35mm from center, 10mm diameter (5mm radius)
        bg_dist_mm = 35
        bg_radius_mm = 5
        # Offset is CCW-positive in math space; image indexing is CW-positive, so subtract for sampling.
        bg_angle_deg = self.ROI_ANGLES[0] - self.angle_offset  # Use first angle for background
        bg_angle_rad = math.radians(bg_angle_deg)

        # Convert background ROI location to pixels
        bg_dist_px = bg_dist_mm / self.pixel_spacing
        bg_radius_px = bg_radius_mm / self.pixel_spacing
        bg_x = cx + bg_dist_px * math.cos(bg_angle_rad)
        bg_y = cy + bg_dist_px * math.sin(bg_angle_rad)

        # Create background mask and extract values
        mask_bg = self._circular_roi_mask(self.image.shape, (bg_x, bg_y), bg_radius_px)
        bg_vals = self.image[mask_bg]

        # Compute background statistics
        if bg_vals.size < 10:
            raise ValueError("Insufficient background pixels for analysis")
        mean_bg = float(np.mean(bg_vals))
        std_bg = float(np.std(bg_vals))

        # Process each ROI
        for roi_name, roi_specs in self.ROI_SETTINGS.items():
            # Extract ROI specifications
            angle_idx = roi_specs["angle_idx"]
            # Offset is CCW-positive in math space; image indexing is CW-positive, so subtract for sampling.
            angle_deg = self.ROI_ANGLES[angle_idx] - self.angle_offset
            radius_mm = roi_specs["radius_mm"]

            # Convert to pixels
            distance_px = self.ROI_DISTANCE_MM / self.pixel_spacing
            radius_px = radius_mm / self.pixel_spacing

            # Calculate ROI center position (match CTP401 convention)
            angle_rad = math.radians(angle_deg)
            x_full = cx + distance_px * math.cos(angle_rad)
            y_full = cy + distance_px * math.sin(angle_rad)

            # Create circular ROI mask
            mask = self._circular_roi_mask(self.image.shape, (x_full, y_full), radius_px)
            vals = self.image[mask]

            # Skip if insufficient data
            if vals.size < 10:
                continue

            # Compute ROI statistics
            mean_signal = float(np.mean(vals))
            std_signal = float(np.std(vals))

            # Contrast: percentage difference from background
            # Positive = brighter than background, Negative = darker
            contrast = ((mean_signal - mean_bg) / mean_bg) * 100.0

            # Contrast-to-Noise Ratio: measures detectability
            # Higher CNR means the ROI is more visible
            cnr = abs(mean_signal - mean_bg) / (std_bg + 1e-8)

            # Store results for this ROI
            # Keep deltas as floats for accurate distance calculation
            x_delta = float(x_full - cx)
            y_delta = float(y_full - cy)
            r_delta = float((x_delta**2 + y_delta**2)**0.5)

            results[f'roi_{roi_name}mm'] = {
                'x': int(x_full),
                'y': int(y_full),
                'r': float(radius_px),
                'x_delta': x_delta,
                'y_delta': y_delta,
                'r_delta': r_delta,
                'angle': float(angle_deg),
                'mean': mean_signal,
                'std': std_signal,
                'bg_mean': mean_bg,
                'bg_std': std_bg,
                'cnr': float(cnr),
                'contrast': float(contrast),
            }

        # Store and return summary
        self.results = {'n_detected': len(results), 'blobs': results}

        if verbose:
            print(f"Low-contrast analysis: {len(results)} ROIs detected")
            for roi_name, roi_data in results.items():
                print(f"  {roi_name}: CNR={roi_data['cnr']:.2f}, Contrast={roi_data['contrast']:.1f}%")
            print()

        return self.results

    def to_dict(self) -> Dict[str, Any]:
        """
        Return JSON-compatible results dictionary.

        Returns:
            Dictionary with all analysis results.
        """
        return self.results

    def get_results_summary(self) -> Dict[str, str]:
        """
        Get a formatted summary of analysis results.

        Returns:
            Dictionary with key measurements formatted as strings.
        """
        if not self.results:
            raise ValueError("Analysis must be run before getting results")

        summary = {}
        summary["ROIs Detected"] = str(self.results['n_detected'])

        # Add CNR for each ROI
        for roi_name, roi_data in self.results['blobs'].items():
            diameter = roi_name.replace('roi_', '').replace('mm', '')
            summary[f"{diameter}mm CNR"] = f"{roi_data['cnr']:.2f}"

        return summary

    def get_plot_data(self) -> Dict[str, Any]:
        """
        Get data needed for plotting visualizations.
        """
        # (Plotter integration handled elsewhere)
        return self.results
