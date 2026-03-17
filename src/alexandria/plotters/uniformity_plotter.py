"""
Uniformity Plotter

Creates comprehensive visualization plots for uniformity analysis results.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Optional
from ..utils import compute_phantom_boundary


class UniformityPlotter:
    """
    Plotter for UniformityAnalyzer results.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.results = analyzer.analyze()

    def _add_roi_box(self, ax, center_xy, size, label, color="yellow", above=True):
        cx, cy = center_xy
        half   = size / 2
        rect = patches.Rectangle((cx - half, cy - half), size, size, linewidth=1.5, edgecolor=color, facecolor='none')
        ax.add_patch(rect)
        stats = self.results[label.lower()]
        text = f"{stats['mean']:.1f} ± {stats['std']:.1f}"
        if above:
            ax.text(cx, cy - 2 * half, text, color=color, ha="center", va="bottom", fontsize=9, bbox=dict(facecolor="black", alpha=0.4, pad=2))
        else:
            ax.text(cx, cy + 2* half, text, color=color, ha="center", va="top", fontsize=9, bbox=dict(facecolor="black", alpha=0.4, pad=2))

    def plot(self):
        img    = self.analyzer.image
        cx, cy = self.analyzer.center
        fig, axes = plt.subplots(3, 2, figsize=(12, 15))
        ax_img    = axes[0, 0]
        ax_hist   = axes[0, 1]
        ax_bar    = axes[1, 0]
        ax_box    = axes[1, 1]
        ax_prof   = axes[2, 0]
        ax_metric = axes[2, 1]
        ax_img.imshow(img, cmap="gray")
        ax_img.set_title("Uniformity Analysis")
        ax_img.set_axis_off()

        pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
        boundary = getattr(self.analyzer, 'boundary', None)
        if boundary and 'x' in boundary and 'y' in boundary:
            boundary_x = np.array(boundary['x'])
            boundary_y = np.array(boundary['y'])
        else:
            _, (boundary_x, boundary_y) = compute_phantom_boundary(img, self.analyzer.center, pixel_spacing)
        if len(boundary_x) > 0:
            ax_img.plot(boundary_x, boundary_y, 'r-', linewidth=1.5, alpha=0.5)

        if hasattr(self.analyzer, 'roi_offset_mm') and pixel_spacing:
            analysis_radius_px = self.analyzer.roi_offset_mm / pixel_spacing
            t = np.linspace(0, 2*np.pi, 100)
            analysis_x = analysis_radius_px * np.cos(t) + self.analyzer.center[0]
            analysis_y = analysis_radius_px * np.sin(t) + self.analyzer.center[1]
            ax_img.plot(analysis_x, analysis_y, 'c--', linewidth=1.0, alpha=0.4)

        ax_img.plot(cx, cy, 'r+', markersize=15, markeredgewidth=2)
        cx, cy = self.analyzer.center
        size   = self.analyzer.roi_size
        offset = self.analyzer.roi_offset
        centers = {"centre": (cx, cy), "north": (cx, cy - offset), "south": (cx, cy + offset), "east": (cx + offset, cy), "west": (cx - offset, cy)}
        roi_colors = { "centre" : "purple", "north"  : "blue", "south"  : "orange", "east"   : "green", "west"   : "red" }
        legend_handles = []
        labels = list(centers.keys())
        means = []
        sems = []
        roi_datas = []
        for label, coord in centers.items():
            color = roi_colors.get(label, "white")
            if label == "centre" or label == "south":
                self._add_roi_box(ax_img, coord, size, label, color, above=False)
            else:
                self._add_roi_box(ax_img, coord, size, label, color)
            cx_roi, cy_roi = coord
            half = size / 2
            roi_data = img[int(cy_roi - half):int(cy_roi + half), int(cx_roi - half):int(cx_roi + half)].flatten()
            roi_datas.append(roi_data)
            ax_hist.hist(roi_data, histtype='step', color=color, linewidth=3, label=label)
            mean_val = self.results[label.lower()]['mean']
            ax_hist.axvline(mean_val, color=color, linestyle='--', linewidth=2)
            legend_handles.append(plt.Line2D([0], [0], color=color, linewidth=2, label=f'{label} (mean: {mean_val:.1f})'))
            means.append(mean_val)
            roi_attr = getattr(self.analyzer, f'm{label[0]}')
            n        = roi_attr.size
            std_val  = self.results[label.lower()]['std']
            sem      = std_val / np.sqrt(n)
            sems.append(sem)

        ax_hist.set_title("ROI Histograms (Overlaid)")
        ax_hist.set_xlabel('HU')
        ax_hist.set_ylabel('Counts')
        ax_hist.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1, 1))
        ax_hist.grid(True, alpha=0.3)

        table_data = [['ROI', 'Mean (HU)', 'Std (HU)', 'SEM (HU)']]
        for label, mean_val, sem in zip(labels, means, sems):
            std_val = self.results[label.lower()]['std']
            table_data.append([label.title(), f"{mean_val:.1f}", f"{std_val:.1f}", f"{sem:.2f}"])
        table = ax_bar.table(cellText=table_data, cellLoc='center', loc='center', colWidths=[0.25, 0.25, 0.25, 0.25])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        for i in range(len(table_data[0])):
            cell = table[(0, i)]
            cell.set_facecolor('#4CAF50')
            cell.set_text_props(weight='bold', color='white')
        for i in range(1, len(table_data)):
            roi_label = labels[i-1]
            for j in range(len(table_data[0])):
                cell = table[(i, j)]
                if i % 2 == 0:
                    cell.set_facecolor('#f0f0f0')
                if j == 0:
                    cell.set_text_props(color=roi_colors[roi_label], weight='bold')

        ax_bar.set_title('ROI Statistics', fontsize=12, weight='bold', pad=10)
        ax_box.boxplot(roi_datas, labels=labels, patch_artist=True)
        ax_box.set_title('ROI Boxplots')
        ax_box.set_ylabel('HU')
        for patch, color in zip(ax_box.patches, [roi_colors[l] for l in labels]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        central_range = 300
        half_range = central_range//2
        center_y = img.shape[0] // 2
        start_y = max(0, center_y - half_range)
        end_y = min(img.shape[0], center_y + half_range)
        vertical_profile = img[start_y:end_y, int(cx)]
        center_x = img.shape[1] // 2
        start_x = max(0, center_x - half_range)
        end_x = min(img.shape[1], center_x + half_range)
        horizontal_profile = img[int(cy), start_x:end_x]
        ax_prof.plot(vertical_profile, label='Vertical (central {}px)'.format(central_range), color='blue')
        ax_prof.plot(horizontal_profile, label='Horizontal (central {}px)'.format(central_range), color='red')
        vert_center_idx = int(cy) - start_y
        horiz_center_idx = int(cx) - start_x
        ax_prof.axvline(vert_center_idx, color='blue', linestyle='--', linewidth=2, label='Vertical center (y={:.0f})'.format(cy))
        ax_prof.axvline(horiz_center_idx, color='red', linestyle='--', linewidth=2, label='Horizontal center (x={:.0f})'.format(cx))
        ax_prof.set_title('Center Profiles (Central 300 Pixels)')
        ax_prof.set_xlabel('Pixel position (relative)')
        ax_prof.set_ylabel('HU')
        ax_prof.legend()
        ax_prof.grid(True, alpha=0.3)

        uni = self.results["uniformity"]
        ax_metric.text(0.5, 0.5, f"Uniformity: {uni:.2f} %", ha='center', va='center', fontsize=16, transform=ax_metric.transAxes)
        ax_metric.axis('off')

        fig.tight_layout()
        return fig
