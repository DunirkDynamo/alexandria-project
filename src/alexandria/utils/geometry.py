"""
Geometry Utilities for CatPhan Phantom Analysis

Provides functions for finding phantom centers, rotations, and geometric measurements.
"""

import numpy as np
from scipy.interpolate import interpn
from scipy.signal import find_peaks, peak_widths
from typing import Tuple, List, Optional, Union


class CatPhanGeometry:
    """Legacy class-based geometry helpers. Prefer the module-level functions for new code."""

    @staticmethod
    def find_center(image: np.ndarray, threshold: float = 400) -> Tuple[List[float], List[np.ndarray]]:
        """
        Locate the phantom centre from a 2-D CT image using threshold crossings.

        Scans the horizontal and vertical profiles through the image midpoint and
        finds the first/last pixel that exceeds ``threshold``.  If no crossings
        are found at the requested threshold, retries at 300 HU.

        Args:
            image    : 2-D CT image array (HU values).
            threshold: HU threshold used to detect the phantom edge (default 400).

        Returns:
            center  : [x, y] coordinates of the phantom centre in pixels.
            boundary: [outer_x, outer_y] arrays tracing the estimated circular
                      phantom boundary.
        """
        # Image dimensions and midpoint pixel for profile extraction
        sz       = np.array(image.shape)          # [rows, cols]
        matrix_c = (np.round(sz[0]/2), np.round(sz[1]/2))  # approximate image centre pixel

        # 1-D profiles through the centre row and centre column
        px = image[int(matrix_c[0]), :]   # horizontal profile (varies along columns)
        py = image[:, int(matrix_c[1])]   # vertical profile (varies along rows)

        # Shift by 1 pixel so the returned coordinate sits inside the phantom edge
        offset = 1

        # Find first and last pixels exceeding the threshold in each profile.
        # These are the phantom edges; if threshold fails, fall back to 300 HU.
        try:
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset   # left edge (col)
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset   # top edge (row)
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset   # right edge (col)
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset   # bottom edge (row)
        except StopIteration:
            # Primary threshold found no crossings; retry at the lower fallback
            threshold = 300
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset

        # Phantom width and height in pixels derived from the edge positions
        szx = x2 - x1   # horizontal extent (cols)
        szy = y2 - y1   # vertical extent (rows)

        # Centre is the midpoint between the left/right and top/bottom edges
        center = [(x1 + x2) / 2, (y1 + y2) / 2]

        # Approximate phantom radius as the average of the two semi-axes (crude but
        # sufficient for a visual boundary overlay)
        outer_r = (szx + szy) / 4

        # Parametric circle for the boundary overlay
        t       = np.linspace(0, 2*np.pi, 100)   # angle parameter
        outer_x = outer_r * np.cos(t) + center[0]
        outer_y = outer_r * np.sin(t) + center[1]
        return center, [outer_x, outer_y]

    @staticmethod
    def select_optimal_ctp528_slices(dicom_set: List, target_index: int, search_range: int = 2) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Select the optimal CTP528 slice image via 3-slice averaging.

        Loads ``2 * search_range + 1`` slices centred on ``target_index``, traces
        a semicircular line-pair profile through each, identifies the slice with the
        highest mean profile intensity, then returns a pixel-averaged image from the
        three slices surrounding that peak (or two slices at the boundary).

        Args:
            dicom_set   : Ordered list of DICOM dataset objects covering the CTP528 module.
            target_index: Index of the intended CTP528 slice within ``dicom_set``.
            search_range: Number of slices either side of ``target_index`` to consider
                          when selecting the optimal slice (default 2).

        Returns:
            im    : Averaged 2-D image array.
            means : Mean profile intensity for each candidate slice.
            z_mean: Mean slice offset (relative to ``target_index``) of the averaged slices.
        """
        z       = target_index                                    # shorthand for the target slice index
        offsets = list(range(-search_range, search_range + 1))   # relative offsets: e.g. [-2,-1,0,1,2]
        imgs    = [dicom_set[z + o].pixel_array for o in offsets] # raw pixel arrays for each candidate slice
        n       = len(imgs)                                        # total number of candidate slices

        # Image geometry taken from the target slice
        sz    = (dicom_set[z].Rows, dicom_set[z].Columns)  # image matrix size (rows, cols)
        space = dicom_set[z].PixelSpacing                   # in-plane pixel spacing [mm/px]
        c     = (int(sz[0] / 2), int(sz[1] / 2))           # image centre pixel (row, col)

        # Semicircular trace through the CTP528 line-pair region.
        # lp_r is the radius of the trace arc in mm (47 mm ≈ line-pair ring radius).
        lp_r  = 47                                                       # trace radius [mm]
        tfine = np.linspace(0, np.pi, 500)                               # angle samples over a semicircle
        lp_b  = lp_r / space[0] * np.cos(tfine) + c[0]                  # x pixel coords of trace arc
        lp_a  = lp_r / space[1] * np.sin(tfine) + c[1]                  # y pixel coords of trace arc

        # Coordinate grids for interpn (scaled by pixel spacing to give physical coords)
        x = np.linspace(0, (sz[0] - 1) / 2, sz[0])   # row coordinate axis
        y = np.linspace(0, (sz[1] - 1) / 2, sz[1])   # col coordinate axis

        # Sample the line-pair trace profile for every candidate slice
        profiles = []
        for img in imgs:
            f = np.zeros(len(lp_a))   # intensity profile along the trace arc
            for i in range(len(lp_a)):
                f[i] = interpn((x, y), img, [lp_a[i] * space[0], lp_b[i] * space[1]])
            profiles.append(f)

        # Identify the candidate slice with the highest mean trace intensity —
        # the CTP528 module has visible line-pair structure that lifts the mean
        means = np.array([np.mean(f) for f in profiles])  # mean profile intensity per candidate slice
        tmp   = int(np.argmax(means))                      # index of the highest-intensity candidate

        # Build a selection mask for the three slices centred on the best candidate.
        # Handles edge cases where the best slice is at the boundary of the search range.
        idx = np.zeros(n)   # 1 = include this slice in the average, 0 = exclude
        try:
            idx[tmp - 1] = 1
            idx[tmp]     = 1
            idx[tmp + 1] = 1
        except IndexError:
            # Best slice is at the very start or end of the candidate window
            if tmp == 0:
                idx[0] = 1
                idx[1] = 1
            elif tmp == n - 1:
                idx[n - 2] = 1
                idx[n - 1] = 1
            else:
                # Unexpected case — return the single best slice unaveraged
                return imgs[tmp].astype(float), means, 0.0

        # Accumulate the selected slices and compute the pixel-wise average
        im         = np.zeros(sz, dtype=float)   # accumulated image sum
        z_mean_list = []                          # offsets of the included slices (relative to target)
        for i, img in enumerate(imgs):
            if idx[i]:
                im += np.array(img, dtype=float)
                z_mean_list.append(offsets[i])
        im     /= float(np.sum(idx))              # normalise to get the pixel average
        z_mean  = float(np.mean(z_mean_list))     # mean z offset of the averaged slices

        return im, means, z_mean

    @staticmethod
    def calculate_slice_thickness(image: np.ndarray, pixel_spacing: float, center: Tuple[float, float]) -> float:
        """
        Measure slice thickness via FWHM of the wire-ramp profile.

        Extracts a strip ROI offset from ``center`` toward the wire ramp, identifies
        the column with the highest integrated signal, locates the ramp peak with
        :func:`scipy.signal.find_peaks`, and converts the FWHM to mm using the
        23-degree ramp angle.

        Args:
            image        : 2-D CT image array containing the wire-ramp feature.
            pixel_spacing: In-plane pixel size in mm.
            center       : (x, y) phantom centre in pixels.

        Returns:
            Slice thickness in mm.
        """
        # Unpack integer pixel coordinates of the phantom centre
        cx, cy = np.int32(center)   # cx = column, cy = row

        # ROI offsets relative to the phantom centre (pixels).
        # The wire ramp is positioned below and to the right of centre.
        roit =  80   # right boundary of the ROI strip (offset from cx)
        roib =  70   # left boundary of the ROI strip (offset from cx)
        roil = -30   # top boundary of the ROI strip (offset from cy)
        roir =  30   # bottom boundary of the ROI strip (offset from cy)

        # Extract the rectangular strip containing the wire ramp
        profs = image[cy + roil:cy + roir, cx + roib:cx + roit]

        # Find the column within the strip that has the strongest total signal —
        # this is the column most likely aligned with the wire ramp peak
        idx_prof = np.argmax(np.sum(profs, axis=0))   # column index of the peak-signal column

        # Peak detection threshold: midpoint between max and min of the profile
        h = (np.max(profs) + np.min(profs)) / 2   # amplitude threshold for find_peaks

        # Locate the ramp peak along the selected column profile
        peaks, _ = find_peaks(profs[:, idx_prof], height=h)

        # Measure the FWHM of the peak; rel_height=0.5 means 50% of peak height
        peaks_results = peak_widths(profs[:, idx_prof], peaks, rel_height=0.5)

        # Convert FWHM from pixels to mm using the 23-degree wire-ramp angle
        # FWHM_mm = FWHM_px * pixel_spacing * sin(23°)
        fwhm = peaks_results[0] * np.sin(np.deg2rad(23)) * pixel_spacing   # slice thickness [mm]
        return float(fwhm[0])


def circular_roi_mask(shape: Tuple[int, int], center: Tuple[float, float], radius: float) -> np.ndarray:
    """
    Create a boolean mask for a circular ROI.

    Args:
        shape : Image dimensions as (height, width).
        center: (x, y) centre of the circle in pixels.
        radius: Circle radius in pixels.

    Returns:
        Boolean array of ``shape``; True inside the circle, False outside.
    """
    # Unpack image dimensions for the coordinate grid
    ny, nx = shape
    # Build a meshgrid of pixel coordinates — Y varies along rows, X along columns
    Y, X = np.ogrid[:ny, :nx]
    cx, cy = center   # circle centre: cx = column, cy = row
    # Euclidean distance from every pixel to the circle centre
    dist = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    # True inside the circle, False outside
    return dist <= radius


def compute_phantom_boundary(image: np.ndarray, center: Tuple[float, float], pixel_spacing: Optional[float] = None, threshold: float = -900, fallback_threshold: float = -900) -> Tuple[Tuple[float, float], Tuple[np.ndarray, np.ndarray]]:
    """
    Estimate the circular phantom boundary from a CT image.

    Scans the row and column that pass through ``center``, finds the first and
    last pixel exceeding ``threshold``, and derives the phantom radius as the
    average of the horizontal and vertical extents.  Falls back to
    ``fallback_threshold`` if the primary threshold yields no crossings.  As a
    last resort, if ``pixel_spacing`` is provided, assumes a 100 mm radius.

    Args:
        image             : 2-D CT image array.
        center            : (x, y) phantom centre in pixels.
        pixel_spacing     : Pixel size in mm; used only for the last-resort fallback.
        threshold         : Primary HU threshold for edge detection (default -900).
        fallback_threshold: Secondary threshold tried when the primary fails.

    Returns:
        detected_center: Refined (x, y) centre derived from the threshold crossings.
        boundary       : (outer_x, outer_y) arrays tracing the estimated circular boundary.
        On complete failure returns ``((0, 0), (empty, empty))``.
    """
    # Guard against degenerate inputs
    if image is None or center is None:
        return ((0, 0), (np.array([]), np.array([])))

    # Clamp the centre coordinates to valid pixel indices
    center_row = int(round(center[1]))   # row index corresponding to the y coordinate
    center_col = int(round(center[0]))   # col index corresponding to the x coordinate
    center_row = max(0, min(center_row, image.shape[0] - 1))
    center_col = max(0, min(center_col, image.shape[1] - 1))

    # 1-D profiles through the clamped centre pixel
    px = image[center_row, :]   # horizontal profile (varies along columns)
    py = image[:, center_col]   # vertical profile (varies along rows)

    # Shift by 1 pixel so the returned coordinate sits inside the phantom edge
    offset = 1

    # Find first and last pixels exceeding the threshold in each profile
    try:
        x1 = next(x for x, val in enumerate(px) if val > threshold) + offset   # left edge (col)
        y1 = next(x for x, val in enumerate(py) if val > threshold) - offset   # top edge (row)
        x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset   # right edge (col)
        y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset   # bottom edge (row)
    except StopIteration:
        # Primary threshold yielded no crossings; try the fallback threshold
        try:
            threshold = fallback_threshold
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
        except StopIteration:
            # Both thresholds failed; use a fixed 100 mm radius if pixel spacing is known
            if pixel_spacing:
                radius_px = 100.0 / pixel_spacing   # 100 mm phantom radius in pixels
                t       = np.linspace(0, 2*np.pi, 100)
                outer_x = radius_px * np.cos(t) + center[0]
                outer_y = radius_px * np.sin(t) + center[1]
                return (center, (outer_x, outer_y))
            return ((0, 0), (np.array([]), np.array([])))

    # Phantom extents from the detected edges
    szx = x2 - x1   # horizontal extent in pixels
    szy = y2 - y1   # vertical extent in pixels

    # Centre is the midpoint between the left/right and top/bottom edges
    detected_center = ((x1 + x2) / 2, (y1 + y2) / 2)

    # Approximate boundary radius as the average of the two semi-axes
    outer_r = (szx + szy) / 4

    # Parametric circle for the boundary overlay
    t       = np.linspace(0, 2*np.pi, 100)   # angle parameter
    outer_x = outer_r * np.cos(t) + detected_center[0]
    outer_y = outer_r * np.sin(t) + detected_center[1]
    return (detected_center, (outer_x, outer_y))


def draw_boundary(center: Tuple[float, float], diameter_x_px: Optional[float], diameter_y_px: Optional[float], n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the perimeter points of an elliptical phantom boundary.

    Args:
        center       : (x, y) centre of the boundary in pixels.
        diameter_x_px: Full horizontal diameter in pixels.
        diameter_y_px: Full vertical diameter in pixels.
        n_points     : Number of points used to discretise the perimeter (default 100).

    Returns:
        (x_coords, y_coords) arrays of boundary points.  Returns empty arrays
        if either diameter is ``None``.
    """
    # Return empty arrays immediately if either diameter is missing
    if diameter_x_px is None or diameter_y_px is None:
        return np.array([]), np.array([])
    rx = diameter_x_px / 2   # horizontal semi-axis in pixels
    ry = diameter_y_px / 2   # vertical semi-axis in pixels
    t        = np.linspace(0, 2 * np.pi, n_points)   # angle parameter
    x_coords = rx * np.cos(t) + center[0]             # x perimeter coords
    y_coords = ry * np.sin(t) + center[1]             # y perimeter coords
    return x_coords, y_coords


def find_center_edge_detection(img: np.ndarray, threshold: float = -900, fallback_threshold: float = -900, return_diameters: bool = False):
    """
    Find the phantom centre using threshold crossings along the image midlines.

    Extracts the central row and column of ``img`` and locates the first and
    last pixels above ``threshold``.  The centre is the midpoint of those
    crossings.  If the primary threshold fails, ``fallback_threshold`` is tried.
    If both fail, the geometric image centre is returned.

    Args:
        img               : 2-D CT image array.
        threshold         : Primary HU threshold for phantom edge detection (default -900).
        fallback_threshold: Secondary threshold used when the primary fails.
        return_diameters  : If True, also return the detected horizontal and
                            vertical diameters in pixels.

    Returns:
        (center_row, center_col) when ``return_diameters`` is False.
        (center_row, center_col, diameter_y_px, diameter_x_px) when True.
        Diameter values are ``None`` if detection failed entirely.
    """
    # Image dimensions and midpoint pixel
    sz       = np.array(img.shape)                                    # [rows, cols]
    matrix_c = (int(np.round(sz[0]/2)), int(np.round(sz[1]/2)))      # centre pixel (row, col)

    # 1-D profiles through the centre row and centre column
    px = img[matrix_c[0], :]   # horizontal profile (varies along columns)
    py = img[:, matrix_c[1]]   # vertical profile (varies along rows)

    # Shift by 1 pixel so detected coordinates sit inside the phantom edge
    offset = 1

    # Find first and last pixels exceeding the threshold in each profile
    try:
        x1 = next(x for x, val in enumerate(px) if val > threshold) + offset   # left edge (col)
        y1 = next(x for x, val in enumerate(py) if val > threshold) - offset   # top edge (row)
        x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset   # right edge (col)
        y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset   # bottom edge (row)
    except StopIteration:
        # Primary threshold failed; retry with the fallback
        threshold = fallback_threshold
        try:
            x1 = next(x for x, val in enumerate(px) if val > threshold) + offset
            y1 = next(x for x, val in enumerate(py) if val > threshold) - offset
            x2 = next(x for x, val in reversed(list(enumerate(px))) if val > threshold) + offset
            y2 = next(x for x, val in reversed(list(enumerate(py))) if val > threshold) - offset
        except StopIteration:
            # Both thresholds failed — return the geometric image centre
            if return_diameters:
                return matrix_c[0], matrix_c[1], None, None
            return matrix_c[0], matrix_c[1]

    # Centre is the midpoint of the detected left/right and top/bottom edges
    center_col = (x1 + x2) / 2.0   # column coordinate of the phantom centre
    center_row = (y1 + y2) / 2.0   # row coordinate of the phantom centre

    if return_diameters:
        # Physical diameter of the phantom in both axes (pixels)
        diameter_x = float(x2 - x1)   # horizontal diameter [px]
        diameter_y = float(y2 - y1)   # vertical diameter [px]
        return center_row, center_col, diameter_y, diameter_x
    return center_row, center_col


def find_rotation(image: np.ndarray, center: Optional[Tuple[float, float]], pixel_spacing: Union[float, Tuple[float, float], List[float]], insert_radius_mm: float = 58.5, edge_threshold: float = 100.0, center_threshold: float = 30, iterations: int = 5, profile_length: int = 25, granularity: int = 4, interp_kwargs: Optional[dict] = None, initial_angle_deg: float = 0.0):
    """
    Detect phantom rotation by locating the top and bottom air/insert positions.

    Seeds two search points at the 270° (top) and 90° (bottom) positions on the
    insert ring.  Each point is iteratively refined by fitting horizontal and
    vertical profiles through the insert and finding the midpoint between the
    two sharpest edges.  The rotation angle is derived from the final vector
    between the two refined points.

    The refinement loop is aborted and the seed positions are restored if any
    iteration moves a point by more than ``center_threshold`` pixels, which
    guards against edge-detection failures on noisy or low-contrast images.

    Args:
        image            : 2-D CT image array.
        center           : (x, y) phantom centre in pixels.  Detected automatically if None.
        pixel_spacing    : Pixel size in mm.  Accepts a scalar, or the first element
                           of a sequence (e.g. DICOM PixelSpacing).
        insert_radius_mm : Distance from phantom centre to the insert ring in mm
                           (default 58.5 mm).
        edge_threshold   : Minimum absolute derivative height required to accept a
                           peak as an insert edge (default 100.0 HU/px).  Lower this
                           for smooth or low-contrast reconstructions.
        center_threshold : Maximum allowed shift per iteration in pixels before
                           the refinement is abandoned (default 30 px).
        iterations       : Maximum number of refinement iterations (default 5).
        profile_length   : Half-length of each sampling profile in pixels (default 25).
        granularity      : Number of sub-pixel samples per pixel along the profile
                           (default 3).  Higher values improve precision at the cost of
                           more interpolation calls.
        interp_kwargs    : Keyword arguments forwarded to ``scipy.interpolate.interpn``.
                           Defaults to ``{'bounds_error': False, 'fill_value': 0}``.
        initial_angle_deg: Seed rotation angle in degrees (default 0).  Use when
                           you know the phantom is already rotated significantly.

    Returns:
        rotation_from_y: Rotation angle in degrees (positive = CCW).
        ct             : Refined (x, y) position of the top insert (nominally 270°).
        cb             : Refined (x, y) position of the bottom insert (nominally 90°).
    """
    # Use a mutable default for interp_kwargs rather than a mutable default argument
    if interp_kwargs is None:
        interp_kwargs = {'bounds_error': False, 'fill_value': 0}

    # Accept pixel spacing as a scalar or as the first element of a DICOM sequence
    if isinstance(pixel_spacing, (list, tuple, np.ndarray)):
        space = float(pixel_spacing[0])   # [mm/px]
    else:
        space = float(pixel_spacing)      # [mm/px]

    # Auto-detect the phantom centre if one was not supplied
    if center is None:
        try:
            r_row, r_col = find_center_edge_detection(image)
            center = (float(r_col), float(r_row))   # (x, y) = (col, row)
        except Exception:
            raise ValueError("Center must be provided or detectable via edge detection")

    # Image dimensions needed to build the coordinate grids
    h_img, w_img = image.shape[:2]   # (rows, cols)

    # Insert ring radius converted from mm to pixels
    ring_r = insert_radius_mm / space   # [px]

    # Seed positions on the insert ring at 270° (top) and 90° (bottom).
    # Image y increases downward, so 90° points down and 270° points up.
    _p90  = (ring_r * np.cos(np.radians(90  + initial_angle_deg)) + center[0],
             ring_r * np.sin(np.radians(90  + initial_angle_deg)) + center[1])   # bottom insert seed (x, y)
    _p270 = (ring_r * np.cos(np.radians(270 + initial_angle_deg)) + center[0],
             ring_r * np.sin(np.radians(270 + initial_angle_deg)) + center[1])   # top insert seed (x, y)
    ct = _p270   # current estimate of the top insert position
    cb = _p90    # current estimate of the bottom insert position

    # Coordinate grids for interpn — row indices [0..h-1], col indices [0..w-1]
    x = np.linspace(0, h_img - 1, h_img)   # row axis
    y = np.linspace(0, w_img - 1, w_img)   # col axis

    def _find_insert_center(roi_pos):
        """Refine an insert centre estimate using edge detection on local profiles.

        Samples a horizontal profile (varying x at fixed y = roi_pos[1]) and a
        vertical profile (varying y at fixed x = roi_pos[0]) around the current
        estimate.  The sharpest pair of edges in each profile are taken as the
        insert boundaries, and their midpoint is returned as the new estimate.
        Returns ``roi_pos`` unchanged if fewer than two edges are found.
        """
        # Sub-pixel sample positions spanning ±profile_length pixels around the current estimate
        x_horiz = np.linspace(roi_pos[0] - profile_length, roi_pos[0] + profile_length, profile_length * granularity)  # col coords of horizontal profile
        x_vert  = np.linspace(roi_pos[1] - profile_length, roi_pos[1] + profile_length, profile_length * granularity)  # row coords of vertical profile

        # Sample image intensities along the horizontal profile (fixed row = roi_pos[1])
        # and the vertical profile (fixed col = roi_pos[0])
        prof_h = np.zeros(len(x_horiz))   # intensity values along the horizontal profile
        prof_v = np.zeros(len(x_vert))    # intensity values along the vertical profile
        for i in range(len(x_horiz)):
            prof_h[i] = interpn((x, y), image, [roi_pos[1], x_horiz[i]], **interp_kwargs)
        for i in range(len(x_vert)):
            prof_v[i] = interpn((x, y), image, [x_vert[i], roi_pos[0]], **interp_kwargs)

        # Discrete first derivative — peaks in |dh| and |dv| correspond to insert edges
        dh = np.diff(prof_h)   # derivative of the horizontal profile
        dv = np.diff(prof_v)   # derivative of the vertical profile
        peaks_h, _ = find_peaks(np.abs(dh), height=edge_threshold)   # edge indices along horizontal profile
        peaks_v, _ = find_peaks(np.abs(dv), height=edge_threshold)   # edge indices along vertical profile

        if len(peaks_h) >= 2 and len(peaks_v) >= 2:
            # Convert derivative peak indices back to pixel coordinates.
            # A peak at index i in dh lies between x_horiz[i] and x_horiz[i+1],
            # so the edge position in pixel coordinates is x_horiz[0] + (i + 0.5) * step.
            step_h = x_horiz[1] - x_horiz[0]   # pixel spacing between horizontal samples
            step_v = x_vert[1]  - x_vert[0]    # pixel spacing between vertical samples
            # Offset of the insert centre from roi_pos, in pixel coordinates
            mid_h = np.mean(x_horiz[0] + (peaks_h + 0.5) * step_h) - roi_pos[0]   # horizontal offset [px]
            mid_v = np.mean(x_vert[0]  + (peaks_v + 0.5) * step_v) - roi_pos[1]   # vertical offset [px]
            return (roi_pos[0] + mid_h, roi_pos[1] + mid_v)

        # Fewer than two edges found — return the current estimate unchanged
        return roi_pos
    # Store original seed positions so they can be restored if refinement diverges
    ct_orig = ct   # original top insert seed
    cb_orig = cb   # original bottom insert seed
    ct_old  = ct   # top insert from the previous iteration (used for the shift sanity check)
    cb_old  = cb   # bottom insert from the previous iteration

    # Iterative refinement: repeatedly re-centre the insert estimate
    for _ in range(iterations):
        try:
            ct_new = _find_insert_center(ct)   # refined top insert position
            cb_new = _find_insert_center(cb)   # refined bottom insert position
        except Exception:
            # Interpolation error (e.g. point outside image bounds) — restore seeds and stop
            ct, cb = ct_orig, cb_orig
            break

        # Sanity check: if either point jumped more than center_threshold pixels in any
        # direction, the edge detection probably failed; restore seeds and stop
        if (abs(ct_new[0] - ct_old[0]) > center_threshold or
                abs(ct_new[1] - ct_old[1]) > center_threshold or
                abs(cb_new[0] - cb_old[0]) > center_threshold or
                abs(cb_new[1] - cb_old[1]) > center_threshold):
            ct, cb = ct_orig, cb_orig
            break

        # Accept the update and advance to the next iteration
        ct_old, cb_old = ct_new, cb_new
        ct, cb = ct_new, cb_new

    # Vector pointing from the bottom insert to the top insert
    tx = ct[0] - cb[0]   # horizontal component (positive = top is to the right)
    ty = ct[1] - cb[1]   # vertical component (positive = top is below, because y increases downward)

    # arctan2(-ty, tx) gives the bearing of the vector in standard math orientation
    # (y-axis pointing up).  Subtracting 90° converts from "angle from east" to
    # "angle from north" (i.e. rotation away from the vertical axis).
    rotation_angle  = np.degrees(np.arctan2(-ty, tx))   # bearing of the top-bottom vector [deg]
    rotation_from_y = float(rotation_angle) - 90.0       # rotation relative to vertical [deg]
    return rotation_from_y, ct, cb
