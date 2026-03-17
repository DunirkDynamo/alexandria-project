# Alexandria — Usage Examples

This document provides concise examples for the current `alexandria` package: analyzers, plotters, and wrappers.

Key notes:
- The shared package does not decide application-specific slice selection for MTF (CTP528). Provide your own selection policy when using DICOM series.
- Most analyzers accept either a single `image` (numpy 2D array) or a `dicom_set` (list of pydicom Dataset) with an optional `slice_index`.

---

## 1) Single-image example — Uniformity

```python
from imageio import imread
from alexandria.analyzers.uniformity import UniformityAnalyzer
from alexandria.wrappers.uniformity_wrapper import UniformityModuleReporter

image = imread('ctp_uniformity.png')  # 2D numpy array

# Analyzer directly
ua = UniformityAnalyzer(image=image, pixel_spacing=0.5)
results = ua.analyze()
print('Uniformity %:', results.get('uniformity'))

# Wrapper convenience
reporter = UniformityModuleReporter(image=image, pixel_spacing=0.5)
reporter.analyze()
fig = reporter.plot()
reporter.save_plot('uniformity_report.png')
```

## 2) High-contrast (MTF) example

```python
from alexandria.analyzers.high_contrast import HighContrastAnalyzer

# Single-image usage
hc = HighContrastAnalyzer(image=image, pixel_spacing=0.5)
mtf_results = hc.analyze()
print('MTF50 lp/mm:', mtf_results['mtf_points'].get('MTF50'))

# DICOM-series usage (application must choose slice_index)
import pydicom
dicom_set = [pydicom.dcmread(p) for p in series_paths]
# choose slice_index using your app's policy, e.g. visually or by metric
hc = HighContrastAnalyzer(dicom_set=dicom_set, slice_index=10)
hc_res = hc.analyze()
```

## 3) CTP401 / CTP404 / CTP515 examples

```python
from alexandria.analyzers.ctp401 import CTP401Analyzer
from alexandria.analyzers.ctp404 import CTP404Analyzer
from alexandria.analyzers.ctp515 import CTP515Analyzer

ct401 = CTP401Analyzer(image=image, pixel_spacing=0.5)
res401 = ct401.analyze()
print(res401.get('ROIs'))

ct404 = CTP404Analyzer(image=image, pixel_spacing=0.5)
res404 = ct404.analyze()

ct515 = CTP515Analyzer(image=image, pixel_spacing=0.5)
res515 = ct515.analyze()
print('CTP515 detected:', res515.get('n_detected'))
```

## 4) Plotters (direct usage)

```python
from alexandria.plotters.uniformity_plotter import UniformityPlotter

plotter = UniformityPlotter(ua)
fig = plotter.plot(figsize=(10, 8))
fig.savefig('uniformity_direct_plot.png')
```

Plotters accept an analyzer instance and return a `matplotlib.Figure`. Use them when you want control over figure creation.

## 5) DICOM-series skeleton (ordering & selection belongs to your app)

```python
import pydicom
from alexandria.analyzers.high_contrast import HighContrastAnalyzer

dicom_set = [pydicom.dcmread(p) for p in series_paths]

# Application-level slice selection example (simple skeleton)
def pick_slice_by_metric(dicom_set):
    # Implement your own policy: metadata, center-of-mass, edge contrast, etc.
    return len(dicom_set)//2

slice_index = pick_slice_by_metric(dicom_set)
hc = HighContrastAnalyzer(dicom_set=dicom_set, slice_index=slice_index)
hc_res = hc.analyze()
```

## 6) Wrapper helpers (common methods)

Wrappers expose a small, consistent API:
- `analyze()` — run analysis and return results
- `plot()` — return a `matplotlib.Figure`
- `analyze_and_plot()` — convenience that returns `(results, fig)`
- `save_plot(path)` — save the last figure
- `get_summary()` — short human-readable string

Check the wrapper modules for constructor options; examples live in this file and in `docs/source/contents.md`.

## 7) Quick smoke test

Before building docs or packaging, run:

```powershell
python -c "import alexandria; print(alexandria.__version__)"
```

---

If you want, I can now:

- scaffold a minimal Sphinx/`autodoc` site under `docs/source/`, or
- expand these examples into runnable Jupyter notebooks or small unit tests.

Tell me which and I'll proceed.
