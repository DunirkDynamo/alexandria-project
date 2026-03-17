"""
CTP401 Plotter

Creates comprehensive visualization plots for 4-ROI linearity analysis results.
Displays main image with ROI circles, per-ROI histograms, and 2D heatmaps.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Optional
from ..utils import compute_phantom_boundary


class CTP401Plotter:
    """
    Plotter for AnalyzerCTP401 results.
    """

    def __init__(self, analyzer, vmin: float = None, vmax: float = None):
        self.analyzer = analyzer
        self.results = analyzer.results
        self.vmin = vmin
        self.vmax = vmax
        self.materials = ['LDPE', 'Air', 'Teflon', 'Acrylic']
        self.roi_angles = { 'LDPE': 0, 'Air': 90, 'Teflon': 180, 'Acrylic': 270 }
        self.roi_colors = { 'LDPE': "blue", 'Air': "orange", 'Teflon': "green", 'Acrylic': "red" }
        self.roi_display_names = { 'LDPE': 'ROI0_LDPE', 'Air': 'ROI90_Air', 'Teflon': 'ROI180_Teflon', 'Acrylic': 'ROI270_Acrylic' }

    def plot(self):
        image = getattr(self.analyzer, "image", None)
        if image is None:
            raise ValueError("Analyzer has no image to plot.")

        fig = plt.figure(figsize=(12, 8))
        gs = fig.add_gridspec(4, 3, width_ratios=[3, 1, 1], wspace=0.3, hspace=0.6)
        ax_img = fig.add_subplot(gs[:, 0])
        ax_img.imshow(image, cmap="gray", vmin=self.vmin, vmax=self.vmax)
        ax_img.set_title("CTP401 ROIs and HU Values")
        ax_img.axis("off")

        boundary = getattr(self.analyzer, 'boundary', None)
        if boundary and 'x' in boundary and 'y' in boundary:
            boundary_x = np.array(boundary['x'])
            boundary_y = np.array(boundary['y'])
        else:
            pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
            _, (boundary_x, boundary_y) = compute_phantom_boundary(image, self.analyzer.center, pixel_spacing)
        
        if len(boundary_x) > 0:
            ax_img.plot(boundary_x, boundary_y, 'r-', linewidth=1.5, alpha=0.5, label='Phantom Boundary')

        pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
        if hasattr(self.analyzer, 'material_distance') and pixel_spacing:
            analysis_radius_px = self.analyzer.material_distance / pixel_spacing
            t = np.linspace(0, 2*np.pi, 100)
            analysis_x = analysis_radius_px * np.cos(t) + self.analyzer.center[0]
            analysis_y = analysis_radius_px * np.sin(t) + self.analyzer.center[1]
            ax_img.plot(analysis_x, analysis_y, 'c--', linewidth=1.0, alpha=0.4, label='Analysis Region')

        cx, cy = self.analyzer.center
        ax_img.plot(cx, cy, 'r+', markersize=15, markeredgewidth=2, label='Center')

        legend_handles = []
        rois = self.analyzer.results.get("ROIs", {})
        ny, nx = image.shape[:2]

        for i, material in enumerate(self.materials):
            stats = rois.get(material)
            if stats is None:
                display_name = self.roi_display_names[material]
                ax_hist = fig.add_subplot(gs[i, 1])
                ax_hist.text(0.5, 0.5, "no data", ha='center', va='center')
                ax_hist.set_title(display_name, fontsize=9)
                ax_hist.axis('off')
                ax_heat = fig.add_subplot(gs[i, 2])
                ax_heat.text(0.5, 0.5, "no data", ha='center', va='center')
                ax_heat.set_title(f"{display_name} Heatmap", fontsize=9)
                ax_heat.axis('off')
                continue

            angle_deg = self.roi_angles[material]
            color = self.roi_colors.get(material, "white")
            display_name = self.roi_display_names[material]

            cx = stats.get('center_x', self.analyzer.center[0])
            cy = stats.get('center_y', self.analyzer.center[1])
            radius_px = self.analyzer.roi_radius / self.analyzer.pixel_spacing

            circle = patches.Circle((cx, cy), radius=radius_px, edgecolor=color, facecolor="none", linewidth=2, alpha=0.6)
            ax_img.add_patch(circle)

            label = f"{display_name}: {stats['mean']:.1f} ± {stats['std']:.1f}"
            legend_handles.append(patches.Circle((0, 0), radius=radius_px, edgecolor=color, facecolor="none", alpha=0.6, label=label))

            ax_hist = fig.add_subplot(gs[i, 1])
            try:
                mask = self.analyzer._create_circular_mask(ny, nx, center=(cx, cy), radius=radius_px)
                data = image[mask]
            except Exception:
                data = np.array([])

            if data.size > 0:
                counts, bins = np.histogram(data.flatten(), bins=30)
                bin_centers = (bins[:-1] + bins[1:]) / 2.0
                import matplotlib.colors as mcolors
                darker_color = mcolors.to_rgba(color, alpha=1.0)
                darker_color = tuple(max(0, c - 0.3) for c in darker_color[:3]) + (1.0,)
                ax_hist.bar(bin_centers, counts, width=(bins[1] - bins[0]), color=color, alpha=0.75, align='center', edgecolor=darker_color, linewidth=0.5)
                mean_val = float(stats.get('mean', np.mean(data)))
                std_val = float(stats.get('std', np.std(data)))
                median_val = float(np.median(data))
                ypos = counts.max() * 0.9 if counts.size else 1.0
                ax_hist.axvline(mean_val, color='k', linestyle='--', linewidth=1)
                ax_hist.axvline(median_val, color='k', linestyle=':', linewidth=1)
                ax_hist.errorbar([mean_val], [ypos], xerr=[std_val], fmt='none', ecolor='k', capsize=3, linewidth=1)
            else:
                ax_hist.text(0.5, 0.5, "no data", ha='center', va='center')

            ax_hist.grid(axis='y', linestyle='--', alpha=0.7)
            for spine in ax_hist.spines.values():
                spine.set_linewidth(2)
                spine.set_color('black')
            ax_hist.set_title(display_name, fontsize=9)
            if i == 3:
                ax_hist.set_xlabel('HU', fontsize=8)
            ax_hist.set_ylabel('Counts', fontsize=8)
            ax_hist.tick_params(axis='both', which='major', labelsize=8)

            ax_heat = fig.add_subplot(gs[i, 2])
            if data.size > 0:
                roi_img = np.zeros((int(radius_px * 2), int(radius_px * 2)))
                y_indices, x_indices = np.ogrid[:roi_img.shape[0], :roi_img.shape[1]]
                dist = np.sqrt((x_indices - radius_px) ** 2 + (y_indices - radius_px) ** 2)
                circle_mask = dist <= radius_px
                roi_img[circle_mask] = data[:circle_mask.sum()]
                ax_heat.imshow(roi_img, cmap='gray', vmin=self.vmin, vmax=self.vmax)
            else:
                ax_heat.text(0.5, 0.5, "no data", ha='center', va='center')
            ax_heat.set_title(f"{display_name} Heatmap", fontsize=9)
            ax_heat.axis('off')

        if legend_handles:
            ax_img.legend(handles=legend_handles, loc="upper right", fontsize=8, framealpha=0.9)

        from matplotlib.lines import Line2D
        mean_line = Line2D([0], [0], color='k', linestyle='--', linewidth=1, label='Mean')
        median_line = Line2D([0], [0], color='k', linestyle=':', linewidth=1, label='Median')
        fig.legend(handles=[mean_line, median_line], loc='lower center', ncol=2, fontsize=10, framealpha=0.9)

        return fig
