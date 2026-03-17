"""
CTP404 Plotter

Creates visualization plots for 9-ROI contrast/sensitometry analysis results.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Optional
from ..utils import compute_phantom_boundary


class CTP404Plotter:
    """
    Creates visualization plots for CTP404Analyzer results.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.vmin = None
        self.vmax = None
        if not self.analyzer.results:
            analyzer.analyze()

    def plot(self, figsize: Optional[tuple] = None, vmin: Optional[float] = None, vmax: Optional[float] = None, **kwargs) -> plt.Figure:
        if vmin is not None:
            self.vmin = vmin
        if vmax is not None:
            self.vmax = vmax

        results = self.analyzer.results.get('contrast', [])
        n_rois = max(1, len(results))
        if figsize is None:
            figsize = (12, max(8, int(n_rois * 1.6)))

        fig = plt.figure(figsize=figsize, constrained_layout=True)
        gs = fig.add_gridspec(n_rois, 3, width_ratios=[3, 1, 1], wspace=0.3, hspace=0.6)
        ax_img = fig.add_subplot(gs[:, 0])
        img = self.analyzer.image
        ax_img.imshow(img, cmap='gray', vmin=self.vmin, vmax=self.vmax)
        ax_img.set_title('CTP404 Sensitometry - 9 Material ROIs')
        
        boundary = getattr(self.analyzer, 'boundary', None)
        if boundary and 'x' in boundary and 'y' in boundary:
            boundary_x = np.array(boundary['x'])
            boundary_y = np.array(boundary['y'])
        else:
            pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
            _, (boundary_x, boundary_y) = compute_phantom_boundary(img, self.analyzer.center, pixel_spacing)
        
        if len(boundary_x) > 0:
            ax_img.plot(boundary_x, boundary_y, 'r-', linewidth=1.5, alpha=0.5, label='Phantom Boundary')

        pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
        if hasattr(self.analyzer, 'material_distance') and pixel_spacing:
            analysis_radius_px = self.analyzer.material_distance / pixel_spacing
            t = np.linspace(0, 2*np.pi, 100)
            analysis_x = analysis_radius_px * np.cos(t) + self.analyzer.center[0]
            analysis_y = analysis_radius_px * np.sin(t) + self.analyzer.center[1]
            ax_img.plot(analysis_x, analysis_y, 'c--', linewidth=1.0, alpha=0.4, label='Analysis Region')

        if self.analyzer.center:
            cx, cy = self.analyzer.center
            ax_img.plot(cx, cy, 'r+', markersize=15, markeredgewidth=2)

        legend_handles = []
        ny, nx = img.shape[:2]
        for i, roi in enumerate(results):
            material = roi.get('material', f'ROI{i+1}')
            cx = roi.get('center_x', self.analyzer.center[0])
            cy = roi.get('center_y', self.analyzer.center[1])
            radius_mm = roi.get('roi_radius_mm', getattr(self.analyzer, 'roi_radius', None))
            radius_px = radius_mm / self.analyzer.pixel_spacing if (radius_mm is not None and self.analyzer.pixel_spacing) else self.analyzer.roi_radius / self.analyzer.pixel_spacing
            color = plt.cm.tab10(i % 10)

            circle = patches.Circle((cx, cy), radius=radius_px, edgecolor=color, facecolor='none', linewidth=2, alpha=0.6)
            ax_img.add_patch(circle)
            label = f"ROI {roi.get('roi_number', i+1)} {material}: {roi.get('mean_hu', float('nan')):.1f} ± {roi.get('std_hu', float('nan')):.1f}"
            legend_handles.append(patches.Circle((0, 0), radius=radius_px, edgecolor=color, facecolor='none', alpha=0.6, label=label))

        ax_img.axis('off')
        ax_img.set_xlabel('')

        from matplotlib.lines import Line2D
        mean_line = Line2D([0], [0], color='k', linestyle='--', linewidth=1, label='Mean')
        median_line = Line2D([0], [0], color='k', linestyle=':', linewidth=1, label='Median')

        for i, roi in enumerate(results):
            ax_hist = fig.add_subplot(gs[i, 1])
            cx = roi.get('center_x', self.analyzer.center[0])
            cy = roi.get('center_y', self.analyzer.center[1])
            radius_mm = roi.get('roi_radius_mm', getattr(self.analyzer, 'roi_radius', None))
            radius_px = radius_mm / self.analyzer.pixel_spacing if (radius_mm is not None and self.analyzer.pixel_spacing) else self.analyzer.roi_radius / self.analyzer.pixel_spacing
            try:
                mask = self.analyzer._create_circular_mask((cx, cy), radius_px)
                data = img[mask]
            except Exception:
                data = np.array([])

            if data.size > 0:
                counts, bins = np.histogram(data.flatten(), bins=30)
                bin_centers = (bins[:-1] + bins[1:]) / 2.0
                ax_hist.bar(bin_centers, counts, width=(bins[1] - bins[0]), color=plt.cm.tab10(i % 10), alpha=0.75, edgecolor='k', linewidth=0.3)
                mean_val = float(roi.get('mean_hu', np.mean(data)))
                std_val = float(roi.get('std_hu', np.std(data)))
                median_val = float(np.median(data))
                ypos = counts.max() * 0.9 if counts.size else 1.0
                ax_hist.axvline(mean_val, color='k', linestyle='--', linewidth=1)
                ax_hist.axvline(median_val, color='k', linestyle=':', linewidth=1)
                ax_hist.errorbar([mean_val], [ypos], xerr=[std_val], fmt='none', ecolor='k', capsize=3, linewidth=1)
            else:
                ax_hist.text(0.5, 0.5, 'no data', ha='center', va='center')

            ax_hist.set_title(f"ROI {roi.get('roi_number', i+1)}: {roi.get('material','')}", fontsize=9)
            ax_hist.grid(axis='y', linestyle='--', alpha=0.7)
            ax_hist.tick_params(axis='both', which='major', labelsize=8)

            ax_heat = fig.add_subplot(gs[i, 2])
            if data.size > 0:
                roi_img = np.zeros((int(radius_px * 2), int(radius_px * 2)))
                y_indices, x_indices = np.ogrid[:roi_img.shape[0], :roi_img.shape[1]]
                dist = np.sqrt((x_indices - radius_px) ** 2 + (y_indices - radius_px) ** 2)
                circle_mask = dist <= radius_px
                flat = data.flatten()
                count = circle_mask.sum()
                if flat.size >= count:
                    roi_img[circle_mask] = flat[:count]
                else:
                    temp = np.full(count, np.nan)
                    temp[:flat.size] = flat
                    roi_img[circle_mask] = temp
                ax_heat.imshow(roi_img, cmap='gray', vmin=self.vmin, vmax=self.vmax)
            else:
                ax_heat.text(0.5, 0.5, 'no data', ha='center', va='center')
            ax_heat.set_title('Heatmap', fontsize=9)
            ax_heat.axis('off')

        if legend_handles:
            ax_img.legend(handles=legend_handles, loc='upper right', fontsize=8, framealpha=0.9)

        fig.legend(handles=[mean_line, median_line], loc='lower center', ncol=2, fontsize=10, framealpha=0.9)
        
        return fig
