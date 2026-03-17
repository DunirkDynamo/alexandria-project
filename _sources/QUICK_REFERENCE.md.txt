# Alexandria — Quick Reference

This quick reference summarizes the current `alexandria` package and common usage patterns. It is intentionally concise — see `docs/CONTENTS.md` for a more detailed reference.

Installation / smoke test

```bash
# From the repository root
pip install -e .[docs]

# Smoke test
python -c "import alexandria; print(alexandria.__version__)"
```

Core modules

- `alexandria.analyzers`:
  - `UniformityAnalyzer`, `DetailedUniformityAnalyzer`, `HighContrastAnalyzer`, `CTP401Analyzer`, `CTP404Analyzer`, `CTP515Analyzer`
- `alexandria.wrappers` — module reporters (convenience): `UniformityModuleReporter`, `HighContrastModuleReporter`, `CTP401ModuleReporter`, `CTP404ModuleReporter`, `CTP515ModuleReporter`
- `alexandria.plotters` — visualizers for analyzers (return `matplotlib.Figure`)
- `alexandria.utils` — `geometry`, `image_processing` helpers

Recommended imports

```python
from alexandria.analyzers.uniformity import UniformityAnalyzer
from alexandria.wrappers.uniformity_wrapper import UniformityModuleReporter
from alexandria.analyzers.high_contrast import HighContrastAnalyzer
```

Usage patterns

1) Wrapper workflow (recommended for scripts)

```python
reporter = UniformityModuleReporter(image=img, pixel_spacing=0.5)
results = reporter.analyze()
fig = reporter.plot()
reporter.save_plot('uniformity.png')
print(reporter.get_summary())
```

2) Direct analyzer usage (fine-grained control)

```python
an = HighContrastAnalyzer(image=img, pixel_spacing=0.5)
res = an.analyze()
print(res['mtf_points'])
```

3) DICOM series (application chooses slice selection)

```python
import pydicom
dicom_set = [pydicom.dcmread(p) for p in series_paths]
slice_index = choose_slice_index(dicom_set)  # app policy
an = HighContrastAnalyzer(dicom_set=dicom_set, slice_index=slice_index)
res = an.analyze()
```

Wrapper common API

- `analyze()` → results (dict)
- `plot()` → `matplotlib.Figure`
- `analyze_and_plot()` → `(results, fig)`
- `save_plot(path)` → save last figure
- `get_summary()` → short human summary

Utilities

- `alexandria.utils.geometry` — center finding, ROI masks, boundary builders, rotation helpers.
- `alexandria.utils.image_processing` — profile extraction, filtering, edge detection, noise estimation.

Best practices

- Provide explicit `center` where available to stabilize ROI placement.
- Use DICOM-series mode for higher SNR; keep slice-selection logic in application code.
- Detect rotation once (e.g., with `CTP401Analyzer`) and reuse the angle for other analyzers.
- Call `analyze()` before `plot()`.

Docs and building locally

Sphinx sources live in `docs/source`. To build locally:

```bash
python -m pip install -e .[docs]
python -m sphinx -b html docs/source docs/_build/html
```

If any analyzer import raises errors during autodoc, ensure the editable install succeeded and required extras are installed.

Contact / contribution

Open issues or PRs for missing examples, additional analyzers, or improved docs. The project welcomes contributions that expand phantom coverage or add reusable medical-physics utilities.
