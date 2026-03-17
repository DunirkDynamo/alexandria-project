"""
CTP404 Analyzer - Contrast and Sensitometry Module

Unified analyzer for CTP404 sensitometry analysis, combining functionality
from catphan404 and XVI-CatPhan implementations.
"""

import numpy as np
from typing import Tuple, List, Dict, Optional, Any
from scipy import ndimage


class CTP404Analyzer:
    """
    Unified CTP404 sensitometry analyzer.

    This analyzer computes mean and standard deviation values for the
    standard CatPhan CTP404 contrast module which contains nine circular
    ROIs containing different materials. The implementation supports two
    initialization modes:

    - Single-image mode: caller passes a prepared 2D NumPy image array
      via the ``image`` parameter.
    - DICOM-series mode: caller provides a list of pydicom dataset
      objects and a ``slice_index``; the analyzer will form a simple
      3-slice average to improve SNR before analysis.
    """
    
    # Material names for each ROI (9 ROIs)
    MATERIALS = [
        'Delrin', 'none', 'Acrylic', 'Air', 'Polystyrene',
        'LDPE', 'PMP', 'Teflon', 'Air2'
    ]
    
    # ROI angles (degrees) relative to phantom rotation
    ROI_ANGLES = [0, 30, 60, 90, 120, 180, -120, -60, -90]
    
    def __init__(
        self,
        image: Optional[np.ndarray] = None,
        dicom_set: Optional[List[Any]] = None,
        slice_index: Optional[int] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,
        rotation_offset: float = 0.0,
        roi_radius: float = 3.5,
        material_distance: float = 58.5,
    ):
        # Initialize image source and analysis mode. We accept either a
        # pre-built NumPy image (convenient for unit tests and scripted
        # workflows) or a DICOM series + slice index. In DICOM mode we
        # produce a simple 3-slice average to reduce noise prior to ROI
        # measurements.
        if image is not None:
            self.image = image.astype(float)
            self.mode = 'single'
        elif dicom_set is not None and slice_index is not None:
            self.dicom_set = dicom_set
            self.slice_index = slice_index
            self.image = self._prepare_averaged_image()
            self.mode = 'dicom'
        else:
            raise ValueError("Must provide either 'image' or both 'dicom_set' and 'slice_index'")

        # Validate and store center coordinates. Center must be provided
        # by the caller (or detected externally) because ROI placement
        # depends on an accurate location.
        if center is None:
            raise ValueError("Must provide 'center' as (x, y) tuple")
        self.center = center

        # Pixel spacing: prefer explicit parameter, fall back to DICOM
        # metadata when available. PixelSpacing may be a sequence (e.g.
        # pydicom MultiValue) so we extract the first element.
        if pixel_spacing is not None:
            self.pixel_spacing = pixel_spacing
        elif spacing is not None:
            self.pixel_spacing = spacing
        else:
            if dicom_set is not None:
                self.pixel_spacing = float(dicom_set[slice_index].PixelSpacing[0])
            else:
                raise ValueError("Must provide 'pixel_spacing' or 'spacing'")

        # Store remaining configuration options
        self.rotation_offset = rotation_offset
        self.roi_radius = roi_radius
        self.material_distance = material_distance

        # Prepare containers for results and plotting coordinates
        self.results = {}
        self.roi_coordinates = []
        self.slice_thickness = None
        
    def _prepare_averaged_image(self) -> np.ndarray:
        """
        Create 3-slice averaged image for improved SNR.
        
        Returns:
            Averaged image array
        """
        idx = self.slice_index
        # Safely form a 3-slice average. When the selected index is at
        # the series endpoints we duplicate the adjacent slice to keep
        # the averaging scheme consistent and avoid index errors.
        if idx == 0:
            im1 = self.dicom_set[idx].pixel_array
            im2 = self.dicom_set[idx + 1].pixel_array
            im3 = self.dicom_set[idx + 1].pixel_array  # Repeat neighbor slice
        elif idx == len(self.dicom_set) - 1:
            im1 = self.dicom_set[idx].pixel_array
            im2 = self.dicom_set[idx - 1].pixel_array
            im3 = self.dicom_set[idx - 1].pixel_array  # Repeat neighbor slice
        else:
            im1 = self.dicom_set[idx].pixel_array
            im2 = self.dicom_set[idx + 1].pixel_array
            im3 = self.dicom_set[idx - 1].pixel_array

        # Convert to float before averaging to avoid integer truncation
        averaged_image = (im1.astype(float) + im2.astype(float) + im3.astype(float)) / 3.0
        return averaged_image
    
    def _create_circular_mask(self, center: Tuple[float, float], radius: float) -> np.ndarray:
        """
        Create a boolean circular mask.
        
        Args:
            center: (x, y) center coordinates
            radius: Radius in pixels
        
        Returns:
            Boolean mask array
        """
        # Create a distance map from the requested center and threshold
        # it to form a boolean mask. Using np.ogrid is memory-efficient
        # for large images because it avoids building full coordinate
        # grids.
        h, w = self.image.shape
        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y - center[1])**2)
        return dist_from_center <= radius

    def _compute_roi_circle(self, angle_deg: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute ROI circle coordinates at specified angle.

        Args:
            angle_deg: Angle in degrees from center

        Returns:
            Tuple of (x_coords, y_coords) for circle perimeter
        """
        cx, cy = self.center
        # Convert the analyzer's image-oriented ROI angle into a math-
        # based angle by subtracting the CCW-positive rotation offset.
        # This keeps ROI placement consistent regardless of phantom
        # rotation detected elsewhere in the pipeline.
        angle_rad = np.deg2rad(angle_deg - self.rotation_offset)

        # Convert configured distances from millimeters to pixels using
        # the stored pixel spacing.
        dist_px = self.material_distance / self.pixel_spacing

        # Compute ROI center coordinates in pixel units
        roi_x = cx + dist_px * np.cos(angle_rad)
        roi_y = cy + dist_px * np.sin(angle_rad)

        # Build a sampled circle perimeter for plotting overlays.
        radius_px = self.roi_radius / self.pixel_spacing
        theta = np.linspace(0, 2 * np.pi, 100)
        x_circle = roi_x + radius_px * np.cos(theta)
        y_circle = roi_y + radius_px * np.sin(theta)

        return x_circle, y_circle

    def detect_rotation(self, initial_angle_deg: float = 0.0) -> float:
        """
        Detect phantom rotation using insert positions (delegates to shared utility).

        Args:
            initial_angle_deg: Initial rotation guess in degrees (default 0)

        Returns:
            Rotation angle in degrees (sets and returns `self.rotation_offset`)
        """
        # Defer imports to avoid circular imports
        from alexandria.utils import find_rotation, find_center_edge_detection

        # If center somehow missing, attempt edge-based center finding
        if getattr(self, 'center', None) is None:
            center_row, center_col = find_center_edge_detection(
                self.image,
                threshold=400.0,
                fallback_threshold=-900.0
            )
            self.center = (center_col, center_row)

        rotation_angle, top_pt, bottom_pt = find_rotation(
            self.image,
            self.center,
            self.pixel_spacing,
            insert_radius_mm=self.material_distance,
            edge_threshold=100.0,
            center_threshold=30,
            iterations=5,
            profile_length=25,
            granularity=3,
            interp_kwargs={'bounds_error': False, 'fill_value': 0},
            initial_angle_deg=initial_angle_deg,
        )

        # Store detected rotation and points for potential external use
        self.rotation_offset = float(rotation_angle)
        # store top/bottom air ROI centers for external spatial-scaling use
        self.rotation_top_point = top_pt
        self.rotation_bottom_point = bottom_pt
        # Return the rotation and the top/bottom insert center points.
        # Some callers expect the tuple (angle, top_pt, bottom_pt) while
        # others request only the angle; returning the full tuple keeps
        # the method flexible.
        return float(self.rotation_offset), top_pt, bottom_pt
    
    def analyze(self, verbose: bool = False) -> Dict:
        """
        Perform contrast analysis on all 9 ROIs.
        
        Args:
            verbose: Print progress information
        
        Returns:
            Dictionary containing analysis results
        """
        if verbose:
            print("Analyzing CTP404 contrast module...")

        contrast_results = []
        self.roi_coordinates = []

        # Iterate over the configured ROI angles and material labels and
        # compute mean/std for each circular ROI. The ROI angles are
        # defined in image-oriented degrees; we subtract the detected
        # CCW-positive rotation offset to obtain the final placement.
        for i, (angle, material) in enumerate(zip(self.ROI_ANGLES, self.MATERIALS)):
            cx, cy = self.center
            angle_rad = np.deg2rad(angle - self.rotation_offset)
            dist_px = self.material_distance / self.pixel_spacing

            # Compute ROI center in pixel coordinates
            roi_x = cx + dist_px * np.cos(angle_rad)
            roi_y = cy + dist_px * np.sin(angle_rad)

            radius_mm = self.roi_radius
            radius_px = radius_mm / self.pixel_spacing

            # Create boolean mask for the circular ROI and sample image
            mask = self._create_circular_mask((roi_x, roi_y), radius_px)
            roi_values = self.image[mask]

            # Aggregate statistics for this ROI
            mean_hu = float(np.mean(roi_values))
            std_hu = float(np.std(roi_values))

            contrast_results.append({
                'roi_number': i + 1,
                'material': material,
                'angle_deg': angle,
                'roi_radius_mm': float(radius_mm),
                'mean_hu': mean_hu,
                'std_hu': std_hu,
                'center_x': roi_x,
                'center_y': roi_y
            })

            # Store polygon/perimeter coordinates for plotting overlays
            theta = np.linspace(0, 2 * np.pi, 100)
            x_circle = roi_x + radius_px * np.cos(theta)
            y_circle = roi_y + radius_px * np.sin(theta)
            self.roi_coordinates.append((x_circle, y_circle))

            if verbose:
                print(f"  ROI {i+1} ({material:12s}): {mean_hu:7.1f} ± {std_hu:5.1f} HU")
        
        # Calculate Low Contrast Visibility (LCV)
        # LCV is calculated from specific material pairs
        lcv = self._calculate_lcv(contrast_results)
        
        self.results = {
            'contrast': contrast_results,
            'LCV_percent': lcv,
            'rotation_offset': self.rotation_offset,
            'mode': self.mode
        }
        
        if verbose:
            print(f"\nLow Contrast Visibility: {lcv:.2f}%")
        
        return self.results
    
    def _calculate_lcv(self, contrast_results: List[Dict]) -> float:
        """
        Calculate Low Contrast Visibility from material ROIs.
        
        Args:
            contrast_results: List of ROI analysis results
        
        Returns:
            LCV percentage
        """
        # Find relevant materials for LCV calculation. The classic LCV
        # metric compares two low-contrast inserts (commonly Polystyrene
        # vs Acrylic) and expresses the difference as a percentage.
        polystyrene = None
        acrylic = None

        for roi in contrast_results:
            if roi['material'] == 'Polystyrene':
                polystyrene = roi['mean_hu']
            elif roi['material'] == 'Acrylic':
                acrylic = roi['mean_hu']

        if polystyrene is not None and acrylic is not None:
            # LCV = |HU1 - HU2| / |HU_ref| * 100 — use the acrylic
            # value as a practical reference. If the reference happens
            # to be zero we fall back to an absolute difference to avoid
            # division-by-zero.
            if acrylic != 0:
                lcv = abs(polystyrene - acrylic) / abs(acrylic) * 100.0
            else:
                lcv = abs(polystyrene - acrylic)
        else:
            lcv = 0.0

        return float(lcv)
    
    def get_results_summary(self) -> str:
        """
        Get formatted summary of analysis results.
        
        Returns:
            Multi-line string summary
        """
        if not self.results:
            return "No analysis results available. Run analyze() first."
        
        lines = ["CTP404 Contrast Analysis Results", "=" * 40]
        
        for roi in self.results['contrast']:
            lines.append(
                f"ROI {roi['roi_number']:2d} ({roi['material']:12s}): "
                f"{roi['mean_hu']:7.1f} ± {roi['std_hu']:5.1f} HU"
            )
        
        lines.append(f"\nLow Contrast Visibility: {self.results['LCV_percent']:.2f}%")
        lines.append(f"Rotation Offset: {self.results['rotation_offset']:.2f}°")
        
        return "\n".join(lines)
    
    def to_dict(self) -> Dict:
        """
        Export results as dictionary.
        
        Returns:
            Results dictionary
        """
        return self.results.copy()
