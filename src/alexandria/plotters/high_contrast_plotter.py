"""
High Contrast Plotter

Creates visualization plots for MTF/resolution analysis results.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Optional


class HighContrastPlotter:
    """
    Creates visualization plots for HighContrastAnalyzer results.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer
        if not hasattr(analyzer, 'nMTF') or analyzer.nMTF is None:
            analyzer.analyze(verbose=False)

    def plot(self, figsize: tuple = None, vmin: float = None, vmax: float = None, **kwargs) -> plt.Figure:
        n_profiles = len(self.analyzer.profiles) if hasattr(self.analyzer, 'profiles') and self.analyzer.profiles else 0
        if figsize is None:
            n_rows = max(n_profiles, 2)
            figsize = (16, max(10, n_rows * 2))

        from matplotlib.gridspec import GridSpec
        n_to_show = n_profiles if n_profiles > 0 else 4
        n_rows = max(n_to_show, 2)
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(n_rows, 2, figure=fig, width_ratios=[1, 1])

        ax_img = fig.add_subplot(gs[:n_rows//2, 0])
        img = self.analyzer.image
        ax_img.imshow(img, cmap='gray', vmin=vmin, vmax=vmax)
        ax_img.set_title('CTP528 Line Pairs & Sampling')

        boundary = getattr(self.analyzer, 'boundary', None)
        if boundary and 'x' in boundary and 'y' in boundary:
            boundary_x = np.array(boundary['x'])
            boundary_y = np.array(boundary['y'])
            if len(boundary_x) > 0:
                ax_img.plot(boundary_x, boundary_y, 'r-', linewidth=1.5, alpha=0.5)

        if hasattr(self.analyzer, 'lp_r_mm') and hasattr(self.analyzer, 'pixel_spacing') and self.analyzer.pixel_spacing:
            analysis_radius_px = self.analyzer.lp_r_mm / self.analyzer.pixel_spacing
            center = self.analyzer.center
            t = np.linspace(0, 2*np.pi, 100)
            analysis_x = analysis_radius_px * np.cos(t) + center[0]
            analysis_y = analysis_radius_px * np.sin(t) + center[1]
            ax_img.plot(analysis_x, analysis_y, 'c--', linewidth=1.0, alpha=0.4)

        if self.analyzer.center:
            cx, cy = self.analyzer.center
            ax_img.plot(cx, cy, 'r+', markersize=15, markeredgewidth=2)

        if self.analyzer.lpx is not None and self.analyzer.lpy is not None:
            ax_img.plot(self.analyzer.lpx, self.analyzer.lpy, '-r', linewidth=1.5, label='Line Pair Path')
            if hasattr(self.analyzer, 'lp_x') and hasattr(self.analyzer, 'lp_y'):
                for xs, ys in zip(self.analyzer.lp_x, self.analyzer.lp_y):
                    ax_img.plot([xs[0], xs[1]], [ys[0], ys[1]], '-g', linewidth=0.8, alpha=0.6)

        ax_img.axis('off')

        ax_mtf = fig.add_subplot(gs[n_rows//2:, 0])
        if self.analyzer.nMTF is not None and self.analyzer.lp_axis is not None:
            ax_mtf.plot(self.analyzer.lp_axis, self.analyzer.nMTF, 'b-', linewidth=2, alpha=0.7, label='nMTF')
            mtf_thresholds = [0.8, 0.5, 0.3, 0.1]
            colors = ['green', 'orange', 'red', 'purple']
            labels = ['80%', '50%', '30%', '10%']
            mtf_points = self.analyzer.mtf_points
            threshold_keys = ['MTF80', 'MTF50', 'MTF30', 'MTF10']
            for threshold, color, label, key in zip(mtf_thresholds, colors, labels, threshold_keys):
                ax_mtf.axhline(threshold, color=color, linestyle='--', alpha=0.3, linewidth=1)
                if key in mtf_points and not np.isnan(mtf_points[key]):
                    ax_mtf.plot(mtf_points[key], threshold, 'o', color=color, markersize=8, mfc='none', label=key, zorder=5)
            ax_mtf.set_xlabel('Spatial Frequency (lp/mm)')
            ax_mtf.set_ylabel('Normalized MTF')
            ax_mtf.set_title('Aggregated Normalized MTF')
            ax_mtf.grid(True, alpha=0.3)
            ax_mtf.legend()
            ax_mtf.set_ylim([0, 1.1])

        if n_to_show > 0 and hasattr(self.analyzer, 'profiles') and self.analyzer.profiles:
            if n_to_show <= 10:
                colors = plt.cm.tab10(np.linspace(0, 1, 10))
            else:
                colors = plt.cm.tab20(np.linspace(0, 1, 20))
            ax_profiles = []
            for i in range(n_to_show):
                if i == 0:
                    ax_prof = fig.add_subplot(gs[i, 1])
                else:
                    ax_prof = fig.add_subplot(gs[i, 1], sharex=ax_profiles[0])
                profile = self.analyzer.profiles[i]
                ax_prof.plot(range(len(profile)), profile, color=colors[i], linewidth=1.5)
                ax_prof.set_ylabel('HU', fontsize=9)
                ax_prof.grid(True, alpha=0.3)
                ax_prof.text(0.02, 0.98, f'LP{i+1} (Profile {i})', transform=ax_prof.transAxes, fontsize=9, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.7, pad=2), color=colors[i], weight='bold')
                if i == n_to_show - 1:
                    ax_prof.set_xlabel('Sample Index')
                else:
                    ax_prof.tick_params(labelbottom=False)
                ax_profiles.append(ax_prof)

        fig.tight_layout()
        return fig
