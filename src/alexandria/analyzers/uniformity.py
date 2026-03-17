"""
Unified Uniformity Analyzer for CatPhan Phantom Analysis

This module combines the uniformity analysis functionality from both
catphan404 and XVI-CatPhan projects, providing a comprehensive analyzer
for CTP486 uniformity module with both single-image and DICOM-series modes.
"""

# Sphinx/autodoc: module docstrings and inline comments expanded for clarity.

import numpy as np
from typing import Optional, Tuple, List, Dict, Any, Callable


class UniformityAnalyzer:
    """
    Analyzer for CT scanner uniformity using the CTP486 module.

    This class evaluates uniformity by measuring mean and standard deviation
    in five fixed ROIs (center, north, south, east, west) relative to the
    phantom center. Supports both single-image analysis and DICOM series with
    3-slice averaging for improved SNR.
    """

    # Region names (standardized across both projects)
    REGIONS = ['centre', 'north', 'south', 'east', 'west']

    def __init__(
        self,
        image: Optional[np.ndarray] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,  # Alias for pixel_spacing
        dicom_set: Optional[List] = None,
        slice_index: Optional[int] = None,
        roi_box_size: float = 15.0,  # in mm
        roi_offset: float = 50.0,     # in mm
        center_finder: Optional[Callable[..., Tuple]] = None,
        center_finder_kwargs: Optional[Dict[str, Any]] = None,
        center_threshold: float = 400.0,
        center_threshold_fallback: float = -900.0,
    ):
        """
        Initialize the UniformityAnalyzer.
        """
        # Mode detection
        self.dicom_mode = dicom_set is not None and slice_index is not None
        
        # DICOM-series mode attributes
        self.dicom_set = dicom_set
        self.slice_index = slice_index
        self.averaged_image = None
        
        # Single-image mode attributes
        if image is not None:
            self.image = np.array(image, dtype=float)
        else:
            self.image = None
        
        self.center = center
        self.center_finder = center_finder
        self.center_finder_kwargs = center_finder_kwargs or {}
        self.center_threshold = center_threshold
        self.center_threshold_fallback = center_threshold_fallback
        
        # Handle spacing parameter aliases
        if pixel_spacing is not None:
            self.pixel_spacing = float(pixel_spacing)
        elif spacing is not None:
            self.pixel_spacing = float(spacing)
        else:
            self.pixel_spacing = None
        
        # ROI parameters (in mm, converted to pixels later)
        self.roi_box_size_mm = float(roi_box_size)
        self.roi_offset_mm = float(roi_offset)
        self.roi_size = None  # Will be set in pixels
        self.roi_offset = None  # Will be set in pixels
        
        # ROI data storage
        self.mc = None  # Center ROI
        self.mn = None  # North ROI
        self.ms = None  # South ROI
        self.me = None  # East ROI
        self.mw = None  # West ROI
        
        # Results storage
        self.results = []
        self.roi_coordinates = None
        self.uniformity_percent = None

    def prepare_image(self):
        """
        Prepare image for analysis.
        """
        if self.dicom_mode:
            # DICOM-series mode: 3-slice averaging
            idx = self.slice_index
            im1 = self.dicom_set[idx].pixel_array
            im2 = self.dicom_set[idx+1].pixel_array
            im3 = self.dicom_set[idx-1].pixel_array
            
            self.averaged_image = (im1 + im2 + im3) / 3
            self.image = self.averaged_image
            
            # Extract center and pixel spacing from DICOM if not provided
            if self.center is None:
                # Center will need to be computed externally or provided
                raise ValueError("Center must be provided for DICOM mode")
            
            if self.pixel_spacing is None:
                self.pixel_spacing = float(self.dicom_set[self.slice_index].PixelSpacing[0])
            
            return self.averaged_image
        else:
            # Single-image mode
            if self.image is None:
                raise ValueError("Image must be provided in single-image mode")
            return self.image

    def _compute_roi_regions(self):
        """
        Compute ROI regions based on center and spacing.
        """
        # Ensure image is prepared
        if self.image is None:
            self.prepare_image()
        
        # Validate pixel spacing
        if self.pixel_spacing is None:
            raise ValueError("Pixel spacing must be provided")
        
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
            self.center = (center_col, center_row)

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
        
        # Convert ROI dimensions from mm to pixels
        self.roi_size = self.roi_box_size_mm / self.pixel_spacing
        self.roi_offset = self.roi_offset_mm / self.pixel_spacing
        
        # Unpack center: center is (x, y) = (col, row)
        cx, cy = self.center
        
        # Compute ROI bounds
        half_size = int(self.roi_size // 2)
        offset = int(self.roi_offset)
        
        # Center ROI
        self.mc = self.image[
            int(cy) - half_size:int(cy) + half_size,
            int(cx) - half_size:int(cx) + half_size
        ]
        
        # North ROI (above center = smaller row, same column)
        self.mn = self.image[
            int(cy) - offset - half_size:int(cy) - offset + half_size,
            int(cx) - half_size:int(cx) + half_size
        ]
        
        # South ROI (below center = larger row, same column)
        self.ms = self.image[
            int(cy) + offset - half_size:int(cy) + offset + half_size,
            int(cx) - half_size:int(cx) + half_size
        ]
        
        # East ROI (right of center = same row, larger column)
        self.me = self.image[
            int(cy) - half_size:int(cy) + half_size,
            int(cx) + offset - half_size:int(cx) + offset + half_size
        ]
        
        # West ROI (left of center = same row, smaller column)
        self.mw = self.image[
            int(cy) - half_size:int(cy) + half_size,
            int(cx) - offset - half_size:int(cx) - offset + half_size
        ]

    def _create_box_mask(self, sz: Tuple[int, int], cx: float, cy: float, roi_sz: float) -> np.ndarray:
        """
        Create a square box mask (alternative method using masks instead of slicing).
        """
        mask = np.zeros(sz)
        x_start = int(cx - roi_sz/2)
        x_end = int(cx + roi_sz/2)
        y_start = int(cy - roi_sz/2)
        y_end = int(cy + roi_sz/2)
        
        # Ensure indices are within bounds
        x_start = max(0, x_start)
        x_end = min(sz[1], x_end)  # sz[1] is width
        y_start = max(0, y_start)
        y_end = min(sz[0], y_end)  # sz[0] is height
        
        # Use proper numpy indexing: [row, col] = [y, x]
        mask[y_start:y_end, x_start:x_end] = 1
        return mask

    def analyze_uniformity(self) -> Tuple[List, np.ndarray, List]:
        """
        Analyze uniformity by measuring HU values in 5 regions.
        """
        # Ensure ROIs are computed
        if self.mc is None:
            self._compute_roi_regions()
        
        # Prepare data structures
        rois = {
            "centre": self.mc,
            "north": self.mn,
            "south": self.ms,
            "east": self.me,
            "west": self.mw
        }
        
        # Calculate statistics for each region
        results = []
        means = []
        
        for region in self.REGIONS:
            roi_data = rois[region]
            mean_val = float(np.mean(roi_data))
            std_val = float(np.std(roi_data))
            results.append([region, mean_val, std_val])
            means.append(mean_val)
        
        # Calculate uniformity metric
        uniformity = (np.max(means) - np.min(means)) / np.max(means) * 100
        results.append(['Uniformity', uniformity, None])
        
        # Store coordinates for plotting
        cx, cy = self.center
        half_sz = int(self.roi_size / 2)
        offset = int(self.roi_offset)
        
        roi_coords = [
            # Centre
            [int(cx) - half_sz, int(cx) + half_sz, int(cy) - half_sz, int(cy) + half_sz],
            # North
            [int(cx) - half_sz, int(cx) + half_sz, int(cy - offset) - half_sz, int(cy - offset) + half_sz],
            # South
            [int(cx) - half_sz, int(cx) + half_sz, int(cy + offset) - half_sz, int(cy + offset) + half_sz],
            # East
            [int(cx + offset) - half_sz, int(cx + offset) + half_sz, int(cy) - half_sz, int(cy) + half_sz],
            # West
            [int(cx - offset) - half_sz, int(cx - offset) + half_sz, int(cy) - half_sz, int(cy) + half_sz],
        ]
        
        # Create composite mask if needed (for plotting)
        if self.dicom_mode:
            sz = self.image.shape
            roi_sz_pixels = self.roi_size
            
            mc_mask = self._create_box_mask(sz, cx, cy, roi_sz_pixels)
            mn_mask = self._create_box_mask(sz, cx, cy - offset, roi_sz_pixels)
            ms_mask = self._create_box_mask(sz, cx, cy + offset, roi_sz_pixels)
            me_mask = self._create_box_mask(sz, cx + offset, cy, roi_sz_pixels)
            mw_mask = self._create_box_mask(sz, cx - offset, cy, roi_sz_pixels)
            
            m_total = mc_mask + mn_mask + ms_mask + me_mask + mw_mask
        else:
            # Create dummy composite mask for single-image mode
            m_total = np.zeros(self.image.shape)
        
        # Store results
        self.results = results
        self.roi_coordinates = roi_coords
        self.uniformity_percent = uniformity
        
        return results, m_total, roi_coords

    def analyze(self, verbose: bool = True) -> Dict[str, Any]:
        """
        Perform the uniformity analysis on five ROIs.
        """
        # Compute ROI regions if not already done
        if self.mc is None:
            self._compute_roi_regions()
        
        # Compute statistics
        rois = {
            "centre": self.mc,
            "north": self.mn,
            "south": self.ms,
            "east": self.me,
            "west": self.mw
        }
        
        results_dict = {}
        means = []

        for name, roi in rois.items():
            mean_val = float(np.mean(roi))
            std_val = float(np.std(roi))
            results_dict[name] = {"mean": mean_val, "std": std_val}
            means.append(mean_val)

        # Overall uniformity metric
        uniformity = (max(means) - min(means)) / max(means) * 100
        results_dict["uniformity"] = uniformity

        return results_dict
