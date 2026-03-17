"""
CTP515 Plotter

Creates comprehensive visualization plots for low-contrast detectability analysis results.
Displays image with color-coded ROIs, dual-axis CNR/Contrast plots, and statistics table.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Optional
from ..utils import compute_phantom_boundary


class CTP515Plotter:
    """
    Plotter for AnalyzerCTP515 (low-contrast detectability) results.

    Creates a 2x2 layout displaying:
      - Image with color-coded ROI circles and background ROI
      - Dual-axis plot of CNR and Contrast vs. ROI diameter
      - Statistics table showing mean, std, CNR, and contrast for each ROI

    ROIs are color-coded by diameter size with adaptive contrast windowing
    to enhance visibility of low-contrast features.

    Args:
        analyzer (AnalyzerCTP515): Completed analyzer instance with results.
    """

    def __init__(self, analyzer):
        """
        Args:
            analyzer (AnalyzerCTP515): Completed analyzer instance with results.
        """
        self.analyzer = analyzer
        self.results  = analyzer.results
        self.image    = analyzer.image
        self.center   = analyzer.center

    def plot(self):
        """
        Generate visualization of low-contrast ROI analysis.
        
        Layout:
          - Top left: Image with ROI overlays
          - Top right: CNR and Contrast vs. ROI Diameter (dual y-axes)
          - Bottom: Statistics table
        """
        fig = plt.figure(figsize=(16, 10))
        
        # Create subplots
        ax_img = plt.subplot2grid((2, 2), (0, 0))
        ax_plot = plt.subplot2grid((2, 2), (0, 1))
        ax_table = plt.subplot2grid((2, 2), (1, 0), colspan=2)
        
        # Left: Image with ROIs (no cropping)
        display_image = self.image
        
        # Apply adaptive contrast focused on background values
        center_crop_size = 100
        h, w             = display_image.shape
        cy_center        = h // 2
        cx_center        = w // 2
        center_sample    = display_image[
            cy_center - center_crop_size//2 : cy_center + center_crop_size//2,
            cx_center - center_crop_size//2 : cx_center + center_crop_size//2
        ]
        bg_mean = np.mean(center_sample)
        bg_std  = np.std(center_sample)
        
        vmin = bg_mean - 3 * bg_std
        vmax = bg_mean + 3 * bg_std
        
        ax_img.imshow(display_image, cmap='gray', vmin=vmin, vmax=vmax)
        ax_img.set_title(f"CTP515 Low-Contrast ROIs (n={self.results['n_detected']})")
        ax_img.axis('off')
        
        # Draw outer phantom boundary (retrieve from analyzer if available)
        boundary = getattr(self.analyzer, 'boundary', None)
        if boundary and 'x' in boundary and 'y' in boundary:
            boundary_x = np.array(boundary['x'])
            boundary_y = np.array(boundary['y'])
        else:
            pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
            _, (boundary_x, boundary_y) = compute_phantom_boundary(display_image, self.center, pixel_spacing)
        
        if len(boundary_x) > 0:
            ax_img.plot(boundary_x, boundary_y, 'r-', linewidth=1.5, alpha=0.5)
        
        # Draw central analysis region (approximate 65mm radius for low-contrast ROIs)
        pixel_spacing = getattr(self.analyzer, 'pixel_spacing', None)
        if pixel_spacing:
            analysis_radius_px = 65.0 / pixel_spacing
            t = np.linspace(0, 2*np.pi, 100)
            analysis_x = analysis_radius_px * np.cos(t) + self.center[0]
            analysis_y = analysis_radius_px * np.sin(t) + self.center[1]
            ax_img.plot(analysis_x, analysis_y, 'c--', linewidth=1.0, alpha=0.4)
        
        # Plot phantom center
        center_col = self.center[0]
        center_row = self.center[1]
        ax_img.plot(center_col, center_row, 'r+', markersize=12, markeredgewidth=2)
        
        # Define color map for different ROI diameters
        color_map = {
            15: '#00ffff',
            9:  '#00ff00',
            8:  '#ffff00',
            7:  '#ff8800',
            6:  '#ff00ff',
            5:  '#ff0000',
        }
        
        diameters       = []
        cnrs            = []
        contrasts       = []
        legend_handles  = []
        legend_labels   = []

        bg_dist_mm   = 35
        bg_radius_mm = 5
        bg_angle_deg = self.analyzer.ROI_ANGLES[0] - self.analyzer.angle_offset
        bg_angle_rad = np.radians(bg_angle_deg)
        spacing      = self.analyzer.pixel_spacing
        bg_dist_px   = bg_dist_mm / spacing
        bg_radius_px = bg_radius_mm / spacing
        cy, cx = int(self.center[1]), int(self.center[0])
        bg_x   = cx + bg_dist_px * np.cos(bg_angle_rad)
        bg_y   = cy + bg_dist_px * np.sin(bg_angle_rad)
        bg_color = 'red'
        bg_circle = patches.Circle((bg_x, bg_y), bg_radius_px, edgecolor=bg_color, 
                                   facecolor='none', linewidth=1, linestyle='-')
        ax_img.add_patch(bg_circle)

        for roi_name, roi_data in self.results['blobs'].items():
            x        = roi_data['x']
            y        = roi_data['y']
            r        = roi_data['r']
            cnr      = roi_data['cnr']
            contrast = roi_data['contrast']
            diameter = float(roi_name.split('_')[1].replace('mm', ''))
            diameters.append(diameter)
            cnrs.append(cnr)
            contrasts.append(contrast)
            color = color_map.get(int(diameter), 'cyan')
            circle = patches.Circle((x, y), r, edgecolor=color, facecolor='none', linewidth=1)
            ax_img.add_patch(circle)
            if int(diameter) not in [int(label.split('mm')[0]) for label in legend_labels]:
                legend_handles.append(patches.Patch(facecolor='none', edgecolor=color, linewidth=1))
                legend_labels.append(f'{int(diameter)}mm')

        legend_handles.append(patches.Patch(facecolor='none', edgecolor=bg_color, linewidth=3, linestyle='-'))
        legend_labels.append('Background')
        if legend_handles:
            ax_img.legend(legend_handles, legend_labels, loc='upper right', fontsize=10, framealpha=0.8, title='ROI Diameter')

        sorted_indices   = np.argsort(diameters)
        diameters_sorted = [diameters[i] for i in sorted_indices]
        cnrs_sorted      = [cnrs[i] for i in sorted_indices]
        contrasts_sorted = [contrasts[i] for i in sorted_indices]

        ax_plot.plot(diameters_sorted, cnrs_sorted, 'bo-', linewidth=3, markersize=8, label='CNR', alpha=0.5)
        ax_plot.set_xlabel('ROI Diameter (mm)', fontsize=12)
        ax_plot.set_ylabel('CNR', fontsize=12, color='blue')
        ax_plot.tick_params(axis='y', labelcolor='blue')
        ax_plot.grid(True, alpha=0.3)
        ax_contrast = ax_plot.twinx()
        ax_contrast.plot(diameters_sorted, contrasts_sorted, 'ro-', linewidth=2, markersize=8, label='Contrast (%)', alpha=0.5)
        ax_contrast.set_ylabel('Contrast (%)', fontsize=12, color='red')
        ax_contrast.tick_params(axis='y', labelcolor='red')
        ax_plot.set_title('CNR and Contrast vs. ROI Diameter')

        lines1, labels1 = ax_plot.get_legend_handles_labels()
        lines2, labels2 = ax_contrast.get_legend_handles_labels()
        ax_plot.legend(lines1 + lines2, labels1 + labels2, loc='lower right', bbox_to_anchor=(1, 0))

        ax_table.axis('off')
        table_data = [['Diameter (mm)', 'Mean (HU)', 'Std (HU)', 'CNR', 'Contrast (%)']]
        roi_list = [(float(name.split('_')[1].replace('mm', '')), name, data) for name, data in self.results['blobs'].items()]
        roi_list.sort(key=lambda x: x[0])
        for diameter, roi_name, roi_data in roi_list:
            table_data.append([
                f"{diameter:.0f}",
                f"{roi_data['mean']:.1f}",
                f"{roi_data['std']:.1f}",
                f"{roi_data['cnr']:.2f}",
                f"{roi_data['contrast']:.2f}"
            ])
        if roi_list:
            bg_data = roi_list[0][2]
            table_data.append([ 'Background', f"{bg_data['bg_mean']:.1f}", f"{bg_data['bg_std']:.1f}", '—', '—' ])

        table = ax_table.table(cellText=table_data, cellLoc='center', loc='center', colWidths=[0.15]*5)
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        for i in range(len(table_data[0])):
            cell = table[(0, i)]
            cell.set_facecolor('#4CAF50')
            cell.set_text_props(weight='bold', color='white')
        for i in range(1, len(table_data)):
            for j in range(len(table_data[0])):
                cell = table[(i, j)]
                if i % 2 == 0:
                    cell.set_facecolor('#f0f0f0')

        ax_table.set_title('ROI Statistics', fontsize=12, weight='bold', pad=10)
        fig.tight_layout()
        return fig
