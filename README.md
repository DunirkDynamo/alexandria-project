# Alexandria - Unified CatPhan Phantom Analysis Library

**Alexandria** is a comprehensive, unified library for analyzing CatPhan CT phantom DICOM images. It combines the best features from multiple CatPhan analysis implementations into a single, well-tested, and maintainable codebase.

## Features

### Comprehensive Analysis Modules

- **UniformityAnalyzer** (CTP486) - Measures image uniformity across 5 ROIs
- **HighContrastAnalyzer** (CTP528) - Line pair resolution and MTF analysis  
- **CTP401Analyzer** (CTP401/404) - Material contrast, HU accuracy, rotation detection
- **CTP515Analyzer** (CTP515) - Low-contrast detectability and CNR measurements

### Dual Operation Modes

Each analyzer supports two modes:
1. **Single-image mode**: Direct analysis of individual CT slices
2. **DICOM-series mode**: Automatic 3-slice averaging for improved SNR

### Advanced Features

- **Automatic rotation detection** using air ROI positions
- **Intelligent slice selection** for optimal image quality
- **Rotation propagation** across analysis modules
- **Flexible initialization** for different workflow patterns
- **JSON-compatible results** for easy integration

## Installation

```bash
# Install from source
cd alexandria
pip install -e .

# Or install from PyPI (when published)
pip install alexandria
```

## Quick Start

### Single-Image Mode

```python
from alexandria import UniformityAnalyzer, HighContrastAnalyzer
import numpy as np

# Load your CT image
image = ...  # 2D numpy array
center = (256, 256)  # (x, y) phantom center
pixel_spacing = 0.488  # mm

# Analyze uniformity
uniformity = UniformityAnalyzer(
    image=image,
    center=center,
    pixel_spacing=pixel_spacing
)
results = uniformity.analyze()
print(f"Uniformity: {results['uniformity']:.2f}%")

# Analyze resolution
high_contrast = HighContrastAnalyzer(
    image=image,
    center=center,
    pixel_spacing=pixel_spacing,
    t_offset_deg=2.5  # rotation offset
)
mtf_results = high_contrast.analyze()
print(f"MTF50: {mtf_results['mtf_50']:.3f} lp/mm")
```

### DICOM-Series Mode

```python
from alexandria import CTP401Analyzer
import pydicom

# Load DICOM series
dicom_files = [...]  # list of DICOM datasets
center = (256, 256)

# Analyze with automatic rotation detection
analyzer = CTP401Analyzer(
    dicom_set=dicom_files,
    slice_index=10,
    center=center,
    use_all_rois=True
)
results = analyzer.analyze(detect_rotation=True)
print(f"Rotation: {results['rotation_angle']:.2f}°")
print(f"LCV: {results['LCV_percent']:.2f}%")
    # Note on angle convention: rotation angles follow image coordinates —
    # 90° corresponds to south (bottom of the image) and 270° corresponds to
    # north (top of the image).

    # Clarification:
    # - The note above describes image coordinates: origin at the top-left
    #   and the y-axis increasing downward (image/screen convention). In
    #   that system positive angles are clockwise (90° → south).
    # - `alexandria.utils.find_rotation` returns a mathematical (Cartesian)
    #   angle in degrees where positive is counter-clockwise (CCW). Callers
    #   that need an image-oriented, clockwise-positive angle for placing
    #   ROIs should convert the returned value (for example by negating it)
    #   before use.
```

### Custom Center-Finding Algorithms

If you want to use your own center-finding method, pass a callable via
`center_finder` when creating an analyzer. This is a simple hook: a function
you provide that Alexandria calls to compute the phantom center instead of the
default edge-based method.

The callable should return either:
- `(row, col)`
- `(row, col, diameter_y_px, diameter_x_px)` if you also want the boundary drawn
  from measured diameters.

```python
from alexandria import CTP401Analyzer

def my_center_finder(image, *, min_hu=-850):
    # Replace with your own logic
    row, col = 250, 260
    return row, col

analyzer = CTP401Analyzer(
    image=ct_image,
    pixel_spacing=0.488,
    center_finder=my_center_finder,
    center_finder_kwargs={"min_hu": -850}
)
results = analyzer.analyze()
```

You can also keep the default finder and adjust its primary and fallback thresholds:

```python
analyzer = CTP401Analyzer(
    image=ct_image,
    pixel_spacing=0.488,
    center_threshold=-980,
    center_threshold_fallback=-900
)
```

## Architecture

Alexandria is designed with modularity and reusability in mind:

```
alexandria/
├── analyzers/          # Core analysis modules
│   ├── uniformity.py   # CTP486 uniformity
│   ├── high_contrast.py # CTP528 resolution/MTF
│   ├── ctp401.py       # CTP401/404 contrast
│   └── ctp515.py       # CTP515 low-contrast
├── plotters/           # Visualization modules
└── utils/              # Shared utilities
```

## Integration with Existing Projects

Alexandria is designed to be a drop-in replacement for existing CatPhan analysis code:

**catphan404 project:**
```python
# Old
from catphan404.uniformity import UniformityAnalyzer

# New
from alexandria import UniformityAnalyzer  # Same API!
```

**XVI-CatPhan project:**
```python
# Old  
from catphan_analysis.modules.ctp486 import CTP486Module

# New
from alexandria import UniformityAnalyzer  # Unified interface
```

## Design Principles

1. **Union of Functionalities**: Each unified analyzer contains the superset of features from both source projects
2. **Backward Compatibility**: Maintains API compatibility with existing code where possible
3. **Mode Flexibility**: Supports both single-image and DICOM-series workflows
4. **Clean Separation**: Analysis logic separated from I/O and visualization
5. **Extensive Documentation**: Clear docstrings and examples

## Contributing

Contributions are welcome! Please ensure:
- Code follows PEP 8 style guidelines
- All tests pass
- Documentation is updated
- Backwards compatibility is maintained

## License

MIT License - See LICENSE file for details

## Credits

Alexandria unifies and extends CatPhan analysis implementations from:
- catphan404 project - Modern modular toolkit
- XVI-CatPhan project - Production-validated workflows

## Version History

### 1.0.0 (2026-01-26)
- Initial release
- Unified analyzers from catphan404 and XVI-CatPhan
- Support for both single-image and DICOM-series modes
- Comprehensive documentation and examples
