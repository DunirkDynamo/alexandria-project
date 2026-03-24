"""
CTP401 Analyzer - Linearity Module (4-ROI)

This module handles the CTP401 linearity analysis with 4 material ROIs:
LDPE, Air, Teflon, and Acrylic positioned at 0°, 90°, 180°, and 270°.

Provides comprehensive analysis including:
- Material ROI contrast measurements (HU accuracy)
- Low Contrast Visibility (LCV)
- Spatial scaling/linearity verification
- Automatic rotation detection using air ROI positions
"""

# Sphinx/autodoc: expanded docstrings and inline comments for maintainability.

import numpy as np
from scipy import ndimage
from scipy.interpolate import interpn
from scipy.signal import find_peaks
from typing import Optional, Tuple, List, Dict, Any, Callable


class CTP401Analyzer:
    """
    Analyzer for CTP401 linearity module (4-ROI version).

    This module measures HU values for 4 material inserts (LDPE, Air, Teflon, Acrylic),
    calculates low contrast visibility, verifies spatial scaling, and can automatically 
    detect phantom rotation.

    Supports both single-image mode and DICOM-series mode with 3-slice averaging.

    Key Features:
    - 4 Material ROI analysis (LDPE at 0°, Air at 90° (south/bottom), Teflon at 180°, Acrylic at 270° (north/top))
    - Automatic rotation detection using air ROI position
    - Low Contrast Visibility (LCV) calculation
    - Spatial scaling verification (X and Y axes)

    Attributes:
        image (np.ndarray): 2D CT image of the module.
        center (tuple): (x, y) coordinates of phantom center in pixels.
        pixel_spacing (float): Pixel spacing in mm.
        results (dict): Analysis results.
    """

    # Material names and angles for 4-ROI mode
    ROI_CONFIG = {
        'LDPE': 0,        # 0 degrees (east/right)
        'Air': 90,        # 90 degrees (south/bottom)
        'Teflon': 180,    # 180 degrees (west/left)
        'Acrylic': 270    # 270 degrees (north/top)
    }

    def __init__(
        self,
        image: Optional[np.ndarray] = None,
        center: Optional[Tuple[float, float]] = None,
        pixel_spacing: Optional[float] = None,
        spacing: Optional[float] = None,
        dicom_set: Optional[List[Any]] = None,
        slice_index: Optional[int] = None,
        roi_radius: float = 3.5,
        material_distance: float = 58.5,
        edge_threshold: float = 100.0,
        center_finder: Optional[Callable[..., Tuple]] = None,
        center_finder_kwargs: Optional[Dict[str, Any]] = None,
        center_threshold: float = -980,
        center_threshold_fallback: float = -900.0,
    ):
        """
        Initialize the CTP401 analyzer (4-ROI linearity module).

        Supports two initialization modes:
        1. Single-image mode: Provide image, center, and pixel_spacing
        2. DICOM-series mode: Provide dicom_set and slice_index

        Args:
            image: 2D CT image of the module (single-image mode).
            center: (x, y) coordinates of phantom center in pixels.
            pixel_spacing: Pixel spacing in mm.
            spacing: Alias for pixel_spacing.
            dicom_set: List of DICOM dataset objects (DICOM-series mode).
            slice_index: Index of the CTP401 slice in the dataset.
            roi_radius: Radius of ROI circles in mm (default: 3.5mm).
            material_distance: Distance from center to ROI centers in mm (default: 58.5mm).
            edge_threshold: Minimum peak height for edge detection (default: 100.0).
                           Tune this for different scanners/protocols if rotation detection fails.
            center_finder: Optional callable to compute center from image.
                           Should return (row, col) or (row, col, diameter_y_px, diameter_x_px).
            center_finder_kwargs: Optional dict of kwargs to pass to center_finder.
            center_threshold: Threshold used by the default center finder.
            center_threshold_fallback: Fallback threshold if the primary fails.
        """
        # Handle image initialization
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
        
        # Handle center - if None, will be computed later
        self.center = center
        
        # Handle pixel spacing (support both parameter names)
        if pixel_spacing is not None:
            self.pixel_spacing = pixel_spacing
        elif spacing is not None:
            self.pixel_spacing = spacing
        else:
            # Try to extract from DICOM
            if dicom_set is not None:
                self.pixel_spacing = float(dicom_set[slice_index].PixelSpacing[0])
            else:
                raise ValueError("Must provide 'pixel_spacing' or 'spacing'")
        
        self.roi_radius                = roi_radius
        self.material_distance         = material_distance
        self.edge_threshold            = edge_threshold
        self.center_finder             = center_finder
        self.center_finder_kwargs      = center_finder_kwargs or {}
        self.center_threshold          = center_threshold
        self.center_threshold_fallback = center_threshold_fallback
        
        # Results storage
        self.results         = {}
        self.roi_coordinates = []
        self.scale           = None
        
    def _prepare_averaged_image(self) -> np.ndarray:
        """
        Create 3-slice averaged image for improved SNR.
        
        Returns:
            Averaged image array
        """
        idx = self.slice_index

        # Prefer averaging distinct slices. If the series is too short to form a
        # three-slice average, fall back to a two-slice average (current + neighbor)
        # instead of counting any slice twice.
        n = len(self.dicom_set)
        if n == 0:
            raise ValueError("dicom_set is empty")

        # Single-slice series: return the slice as float image
        if n == 1:
            return self.dicom_set[0].pixel_array.astype(float)

        # Edge cases: use a two-slice average (current + neighbor)
        if idx <= 0:
            im0 = self.dicom_set[0].pixel_array.astype(float)
            im1 = self.dicom_set[1].pixel_array.astype(float)
            return (im0 + im1) / 2.0
        if idx >= n - 1:
            imm1 = self.dicom_set[n - 2].pixel_array.astype(float)
            imn = self.dicom_set[n - 1].pixel_array.astype(float)
            return (imm1 + imn) / 2.0

        # Typical case: average previous, current, and next slices
        im_prev = self.dicom_set[idx - 1].pixel_array.astype(float)
        im_curr = self.dicom_set[idx].pixel_array.astype(float)
        im_next = self.dicom_set[idx + 1].pixel_array.astype(float)
        return (im_prev + im_curr + im_next) / 3.0
    
    def _create_circular_mask(
        self,
        h: int,
        w: int,
        center: Optional[Tuple[float, float]] = None,
        radius: Optional[float] = None
    ) -> np.ndarray:
        """
        Create a boolean circular mask.

        Args:
            h: Image height.
            w: Image width.
            center: (x, y) coordinates of circle center.
            radius: Circle radius in pixels.

        Returns:
            Boolean mask array.
        """
        if center is None:
            # Default to geometric image center when none provided.
            center = (int(w / 2), int(h / 2))
        if radius is None:
            # If radius not provided, use the largest circle that fits in image.
            radius = min(center[0], center[1], w - center[0], h - center[1])

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y - center[1])**2)
        return dist_from_center <= radius
    
    def analyze(self, t_offset: float = 0.0, verbose: bool = False) -> Dict:
        """
        Perform ROI analysis on the stored image.

        Args:
            t_offset: Rotational offset for ROIs in degrees
            verbose: Print progress information

        Returns:
            Dictionary containing analysis results with:
            - 'ROIs': Dictionary of ROI results (mean, std for each material)
            - 'LCV_percent': Low contrast visibility percentage
            - 'Scale': Spatial scaling factors (scaleX_cm, scaleY_cm)
        """
        if verbose:
            print("Analyzing CTP401 linearity module...")
        
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
        
        image = self.image
        h, w  = image.shape[:2]
        c0    = self.center
        space = self.pixel_spacing

        # ROI parameters in pixels
        r      = self.roi_radius / space          # ROI radius in pixels
        ring_r = self.material_distance / space  # Distance to ROI centers in pixels

        roi_results          = {}
        self.roi_coordinates = []

        # Analyze each of the 4 ROIs
        for material, angle in self.ROI_CONFIG.items():
            # Calculate ROI center
            # Offset is CCW-positive in math space; image indexing is CW-positive, so subtract for sampling.
            angle_rad = np.radians(angle - t_offset)
            cx        = ring_r * np.cos(angle_rad) + c0[0]
            cy        = ring_r * np.sin(angle_rad) + c0[1]

            # Create circular mask
            mask = self._create_circular_mask(h, w, center=(cx, cy), radius=r)

            # Calculate statistics
            roi_mean = float(np.mean(image[mask]))
            roi_std  = float(np.std(image[mask]))
            
            roi_results[material] = {
                'mean': roi_mean,
                'std': roi_std,
                'angle': angle,
                'center_x': cx,
                'center_y': cy
            }

            # Store ROI coordinates for visualization
            t_viz = np.linspace(0, 2*np.pi, 100)
            roi_x = r * np.cos(t_viz) + cx
            roi_y = r * np.sin(t_viz) + cy
            self.roi_coordinates.append((roi_x, roi_y))
            
            if verbose:
                print(f"  {material:8s} (angle {angle:4.0f}°): {roi_mean:7.1f} ± {roi_std:5.1f} HU")

        # Compute Low Contrast Visibility (LCV)
        # LCV = 3.25 * (σ_air + σ_LDPE) / (μ_air - μ_LDPE)
        denominator = roi_results['Air']['mean'] - roi_results['LDPE']['mean']
        if abs(denominator) < 1e-6:
            print(f"⚠️  Warning: Air and LDPE ROIs have identical means ({roi_results['Air']['mean']:.2f} HU)")
            print("    This suggests incorrect ROI positioning or invalid image data.")
            lcv = 0.0  # Set to 0 rather than causing division by zero
        else:
            lcv = 3.25 * (roi_results['Air']['std'] + roi_results['LDPE']['std']) / denominator

        # Compute spatial scaling factors
        scale = self._compute_spatial_scaling(t_offset)

        # Store results
        self.results = {
            'ROIs': roi_results,
            'LCV_percent': float(lcv),
            'Scale': scale,
            'rotation_offset': t_offset,
            'mode': self.mode
        }
        
        if verbose:
            print(f"\nLow Contrast Visibility: {lcv:.2f}%")
            print(f"Spatial Scaling: X={scale['scaleX_cm']:.2f} cm, Y={scale['scaleY_cm']:.2f} cm")

        return self.results
    
    def _compute_spatial_scaling(self, t_offset: float) -> Dict[str, float]:
        """
        Compute spatial scaling/linearity using edge detection.
        
        Args:
            t_offset: Rotation offset in degrees
        
        Returns:
            Dictionary with scaleX_cm and scaleY_cm
        """
        image  = self.image
        c0     = self.center
        space  = self.pixel_spacing
        ring_r = self.material_distance / space
        
        # Extract horizontal and vertical profiles through center
        px             = image[int(round(c0[1])), :].astype(float)
        py             = image[:, int(round(c0[0]))].astype(float)
        profile_length = 26  # pixels to sample around each ROI

        # Calculate positions of ROIs along axes
        # Offset is CCW-positive in math space; image indexing is CW-positive, so subtract for sampling.
        idx_x1 = int(round(ring_r * np.cos(np.radians(0 - t_offset)) + c0[0]))
        idx_x2 = int(round(ring_r * np.cos(np.radians(180 - t_offset)) + c0[0]))
        idx_y1 = int(round(ring_r * np.sin(np.radians(90 - t_offset)) + c0[1]))
        idx_y2 = int(round(ring_r * np.sin(np.radians(270 - t_offset)) + c0[1]))

        # Extract profiles around ROIs
        px1 = px[max(0, idx_x1-profile_length):min(len(px), idx_x1+profile_length)]
        px2 = px[max(0, idx_x2-profile_length):min(len(px), idx_x2+profile_length)]
        py1 = py[max(0, idx_y1-profile_length):min(len(py), idx_y1+profile_length)]
        py2 = py[max(0, idx_y2-profile_length):min(len(py), idx_y2+profile_length)]

        # Take derivatives to find edges
        dpx1, dpx2 = np.diff(px1), np.diff(px2)
        dpy1, dpy2 = np.diff(py1), np.diff(py2)

        # Find edge positions (min and max of derivative)
        minx1, maxx1 = np.argmin(dpx1), np.argmax(dpx1)
        minx2, maxx2 = np.argmin(dpx2), np.argmax(dpx2)
        miny1, maxy1 = np.argmin(dpy1), np.argmax(dpy1)
        miny2, maxy2 = np.argmin(dpy2), np.argmax(dpy2)

        # Calculate scaling factors (expected 10cm between opposite ROIs)
        scaleX1 = (abs(idx_x2 - idx_x1) * space - abs(minx1 - minx2) * space) / 10
        scaleX2 = (abs(idx_x2 - idx_x1) * space - abs(maxx1 - maxx2) * space) / 10
        scaleY1 = (abs(idx_y1 - idx_y2) * space - abs(miny1 - miny2) * space) / 10
        scaleY2 = (abs(idx_y1 - idx_y2) * space - abs(maxy1 - maxy2) * space) / 10

        scaleX = float(np.mean([scaleX1, scaleX2]))
        scaleY = float(np.mean([scaleY1, scaleY2]))

        self.scale = {'scaleX_cm': scaleX, 'scaleY_cm': scaleY}
        return self.scale
    
    def detect_rotation(self, initial_angle_deg: float = 0.0) -> float:
        """
        Detect phantom rotation using material insert positions.
        
        Finds the Air insert (90°, top) and Acrylic insert (270°, bottom) and 
        calculates rotation based on their offset from vertical alignment. 
        Uses iterative refinement with safety threshold to prevent invalid results.
        """
        # Implementation omitted in the copied excerpt to keep file concise for initial migration.
        # The full implementation remains in the original package and will be migrated next.
        return initial_angle_deg
