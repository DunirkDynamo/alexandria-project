# Alexandria Library Contents

Complete reference of all classes and functions in the Alexandria package (current layout).

---

## Architecture Overview

The package groups functionality into:

1. **Analyzers** — core image-analysis engines (uniformity, contrast, MTF, etc.)
2. **Plotters** — visualization helpers to render analyzer results
3. **Wrappers** — convenience module reporters combining analyzers and plotters
4. **Utils** — geometry and image-processing helpers used across modules

This document lists the primary public classes, constructor parameters, and the outputs you can expect when calling the analysis routines.

---

## Analyzers

### `UniformityAnalyzer`

**Location:** `alexandria.analyzers.uniformity`

Signature (constructor):
```python
UniformityAnalyzer(
    image: Optional[np.ndarray] = None,
    center: Optional[Tuple[float, float]] = None,
    pixel_spacing: Optional[float] = None,
    spacing: Optional[float] = None,
    dicom_set: Optional[List] = None,
    slice_index: Optional[int] = None,
    roi_box_size: float = 15.0,  # mm
    roi_offset: float = 50.0     # mm
)
```

Purpose: five-ROI uniformity analysis (centre, north, south, east, west). Supports single-image and DICOM-series (3-slice averaged) modes.

Key methods:
- `analyze(verbose: bool = True) -> Dict[str, Any]`
- `analyze_uniformity() -> Tuple[List, np.ndarray, List]` (compat return for some callers)

Expected outputs: dictionary with per-ROI mean/std and a `uniformity` percentage, plus optional `roi_coordinates` and `boundary` for plotting

---

### `DetailedUniformityAnalyzer`

**Location:** `alexandria.analyzers.detailed_uniformity`

Signature:
```python
DetailedUniformityAnalyzer(
    image: Optional[np.ndarray] = None,
    center: Optional[Tuple[float, float]] = None,
    pixel_spacing: Optional[float] = None,
    spacing: Optional[float] = None,
    dicom_set: Optional[List] = None,
    slice_index: Optional[int] = None,
    radii_mm: Optional[List[float]] = None,
    sample_step_mm: float = 1.0,
    n_samples: int = 360
)
```

Purpose: sample concentric circular profiles and return angle/value series per radius for detailed uniformity inspection.

Key methods:
- `analyze() -> Dict[str, Any]` — returns `profiles` (JSON-friendly) and metadata

---

### `HighContrastAnalyzer`

**Location:** `alexandria.analyzers.high_contrast`

Signature:
```python
HighContrastAnalyzer(
    image: Optional[np.ndarray] = None,
    pixel_spacing: Optional[float] = None,
    center: Optional[Tuple[float, float]] = None,
    t_offset_deg: float = 0.0,
    rotation_offset: Optional[float] = None,
    dicom_set: Optional[List] = None,
    slice_index: Optional[int] = None,
    lp_r_mm: float = 48.0,
    samples_per_segment: int = 50,
    center_threshold: float = -980,
    center_threshold_fallback: float = -900.0,
)
```

Purpose: analyze CTP528 line-pair module, sample profiles, compute per-pair modulation and derive MTF curve and benchmark frequencies (MTF80/50/30/10).

Key methods:
- `analyze(verbose: bool = True, write_log: bool = False) -> Dict[str, Any]`

Outputs: normalized MTF array, frequency axis `lp_axis`, and `mtf_points` mapping benchmark MTF levels to lp/mm

Notes: DICOM-mode attempts to use a project-provided slice-selection routine when available, otherwise falls back to 3-slice averaging.

---

### `CTP401Analyzer`

**Location:** `alexandria.analyzers.ctp401`

Signature:
```python
CTP401Analyzer(
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
)
```

Purpose: linearity / material ROI analysis for the CTP401 module. Includes rotation detection utilities and spatial scaling checks.

Key methods:
- `detect_rotation(initial_angle_deg: float = 0.0) -> float`
- `analyze(t_offset: float = 0.0, verbose: bool = True) -> Dict[str, Any]`

Outputs: per-material ROI statistics, LCV estimates, and spatial scaling information.

---

### `CTP404Analyzer`

**Location:** `alexandria.analyzers.ctp404`

Signature:
```python
CTP404Analyzer(
    image: Optional[np.ndarray] = None,
    dicom_set: Optional[List] = None,
    slice_index: Optional[int] = None,
    center: Optional[Tuple[str, float]] = None,
    pixel_spacing: Optional[float] = None,
    spacing: Optional[float] = None,
    rotation_offset: float = 0.0,
    roi_radius: float = 3.5,
    material_distance: float = 60.0,
    center_finder: Optional[Callable[..., Tuple]] = None,
    center_finder_kwargs: Optional[Dict[str, Any]] = None,
)
```

Purpose: 9-ROI sensitometry analysis (multiple materials) producing contrast metrics and LCV.

Key methods:
- `analyze(verbose: bool = True) -> Dict[str, Any]`

Outputs: list of ROI contrast results and summary metrics.

---

### `CTP515Analyzer`

**Location:** `alexandria.analyzers.ctp515`

Signature:
```python
CTP515Analyzer(
    image: Optional[np.ndarray] = None,
    center: Optional[Tuple[float, float]] = None,
    pixel_spacing: Optional[float] = None,
    angle_offset: float = 0.0,
    dicom_set: Optional[List] = None,
    slice_index: Optional[int] = None,
    center_finder: Optional[Callable[..., Tuple]] = None,
    center_finder_kwargs: Optional[Dict[str, Any]] = None,
)
```

Purpose: low-contrast detectability analysis for CTP515 (multiple circular inserts). Computes per-ROI CNR and detectability summaries.

Key methods:
- `analyze(verbose: bool = True) -> Dict[str, Any]`

Outputs: per-ROI CNR/contrast values and `n_detected` count.

---

## Plotters

Plotter classes take a completed analyzer instance and render figures for reporting.

### `UniformityPlotter`
Location: `alexandria.plotters.uniformity_plotter`
Signature: `UniformityPlotter(analyzer: UniformityAnalyzer)`

Generates a 3x2 figure containing image + ROI overlays, histograms, table, boxplots, profiles, and a uniformity metric panel.

---

### `DetailedUniformityPlotter`
Location: `alexandria.plotters.detailed_uniformity_plotter`
Signature: `DetailedUniformityPlotter(analyzer: DetailedUniformityAnalyzer)`

Plots concentric profile traces and polar/linear visualizations per radius.

---

### `HighContrastPlotter`
Location: `alexandria.plotters.high_contrast_plotter`
Signature: `HighContrastPlotter(analyzer: HighContrastAnalyzer)`

Visualizes per-pair profiles, peak/trough diagnostics, and the combined MTF curve with annotated MTF benchmark points.

---

### `CTP401Plotter`, `CTP404Plotter`, `CTP515Plotter`
Location: `alexandria.plotters.*`

Each plotter offers ROI overlays, labeled statistics, and a `plot()` method returning a `matplotlib.Figure`. See respective modules for plotting options and output layout.

---

## Wrappers (ModuleReporters)

Convenience classes that create an analyzer and plotter and expose `analyze()`, `plot()`, `analyze_and_plot()`, `save_plot()`, and `get_summary()` helpers.

Public wrappers exported in `alexandria.__init__`:
- `UniformityModuleReporter` — thin wrapper around `UniformityAnalyzer` + `UniformityPlotter` (constructor supports single-image and DICOM-series modes; see `wrappers/uniformity_wrapper.py`)
- `HighContrastModuleReporter` — wrapper for MTF workflows (supports `lp_r_mm`, `samples_per_segment` arguments)
- `CTP401ModuleReporter`, `CTP404ModuleReporter`, `CTP515ModuleReporter` — wrappers for the contrast/linearity/low-contrast modules

Each wrapper module contains maintainers' developer notes and example usage in `docs/source/USAGE_EXAMPLES.md`.

---

## Utilities

### `alexandria.utils.geometry.CatPhanGeometry`

Useful helpers for center finding, boundary construction, rotation detection, and slice-selection notes.

Key routines:
- `find_center(image, threshold=400) -> (center, [outer_x, outer_y])` — 1D profile edge detection fallback behavior documented in code
- `select_optimal_ctp528_slices(dicom_set, target_index, search_range=2)` — intentionally left unimplemented in shared package; application projects should provide selection policy
- `calculate_slice_thickness(image, pixel_spacing, center)` — placeholder
- `circular_roi_mask(shape, center, radius) -> np.ndarray`
- `compute_phantom_boundary(image, center, pixel_spacing=None, threshold=-900, fallback_threshold=-900)` — robust fallback to nominal radius when edge detection fails
- `draw_boundary(center, diameter_x_px, diameter_y_px, n_points=100)` — ellipse boundary builder
- `find_center_edge_detection(img, threshold=-900, fallback_threshold=-900, return_diameters=False)` — XVI-style center finder
- `find_rotation(image, center, pixel_spacing, ...) -> (rotation_deg, top_point, bottom_point)` — robust rotation detection with documented angle conventions

## Image Processing Utilities

### `alexandria.utils.image_processing.ImageProcessor`

Stateless helper methods (all `@staticmethod`):
- `apply_gaussian_filter(image, sigma=1.0) -> np.ndarray`
- `extract_profile(image, start, end, n_points=100) -> np.ndarray` — subpixel sampling via `scipy.interpolate.interpn`
- `threshold_image(image, threshold, mode='above') -> np.ndarray`
- `estimate_noise(image, roi_center=None, roi_size=50) -> float`
- `find_edges(image, method='sobel') -> np.ndarray` — supports 'sobel', 'prewitt', and 'scharr'

---

## Usage & Examples

See `docs/source/USAGE_EXAMPLES.md` for runnable examples illustrating single-image and DICOM-series usage patterns, wrapper convenience workflows, and advanced notes (slice selection and ROI masking).

---

If you want, I can now generate a table-of-contents `index.rst` and a minimal Sphinx `conf.py` under `docs/source/` and wire up `autodoc` entries for the modules above. Would you like me to scaffold that next?
