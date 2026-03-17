"""
Detailed Uniformity Plotter

Plots concentric profile values (angle vs HU) and overlayed histograms.
"""

from typing import List, Dict, Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter


class DetailedUniformityPlotter:
    """
    Plotter for DetailedUniformityAnalyzer results.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer
        if not getattr(self.analyzer, "results", None):
            self.analyzer.analyze()

    def _get_profiles(self) -> List[Dict[str, Any]]:
        plot_radii = set(round(float(r), 3) for r in getattr(self.analyzer, "radii_mm", self.analyzer.results.get("plot_radii_mm", [])))
        if getattr(self.analyzer, "profile_data", None):
            return [profile for profile in self.analyzer.profile_data if round(float(profile.get("radius_mm", 0.0)), 3) in plot_radii]
        profiles = []
        for profile in self.analyzer.results.get("profiles", []):
            if round(float(profile.get("radius_mm", 0.0)), 3) not in plot_radii:
                continue
            profiles.append({"radius_mm": profile["radius_mm"], "angles_deg": np.array(profile["angles_deg"], dtype=float), "values": np.array(profile["values"], dtype=float)})
        return profiles

    def _get_all_profiles(self) -> List[Dict[str, Any]]:
        if getattr(self.analyzer, "profile_data", None):
            return self.analyzer.profile_data
        profiles = []
        for profile in self.analyzer.results.get("profiles", []):
            profiles.append({"radius_mm": profile["radius_mm"], "angles_deg": np.array(profile["angles_deg"], dtype=float), "values": np.array(profile["values"], dtype=float)})
        return profiles

    def plot(self, bins: int = 25, figsize: tuple = (12, 18)) -> plt.Figure:
        profiles = self._get_profiles()
        if not profiles:
            raise ValueError("No profile data available for plotting")

        fig, (ax_profile, ax_smooth, ax_resid, ax_hist, ax_mean, ax_img) = plt.subplots(6, 1, figsize=figsize)
        colors = plt.cm.tab10(np.linspace(0, 1, max(3, len(profiles))))

        mean_points = []
        std_points = []
        radius_points = []

        for idx, profile in enumerate(profiles):
            radius = profile["radius_mm"]
            angles = profile["angles_deg"]
            values = profile["values"]
            color = colors[idx % len(colors)]
            ax_profile.plot(angles, values, color=color, linewidth=1.2, label=f"r={radius:.1f}mm")
            window_length = min(17, len(values) - 1 if len(values) % 2 == 0 else len(values))
            if window_length < 5:
                smoothed = values
            else:
                if window_length % 2 == 0:
                    window_length -= 1
                smoothed = savgol_filter(values, window_length=window_length, polyorder=3, mode='interp')
            ax_smooth.plot(angles, smoothed, color=color, linewidth=1.2, label=f"r={radius:.1f}mm")
            mean_val = float(np.mean(values))
            residuals = values - mean_val
            ax_resid.plot(angles, residuals, color=color, linewidth=1.0, label=f"r={radius:.1f}mm")
            ax_hist.hist(values, bins=bins, histtype='step', color=color, linewidth=2, label=f"r={radius:.1f}mm")
            std_val = float(np.std(values))
            if std_val > 0.0:
                x_min = float(np.min(values))
                x_max = float(np.max(values))
                if x_max > x_min:
                    x_vals = np.linspace(x_min, x_max, 200)
                    bin_width = (x_max - x_min) / float(bins)
                    pdf = (1.0 / (std_val * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x_vals - float(np.mean(values))) / std_val) ** 2)
                    ax_hist.plot(x_vals, pdf * len(values) * bin_width, color=color, linestyle='--', linewidth=1.0)
            radius_points.append(radius)
            mean_points.append(mean_val)
            std_points.append(std_val)

        ax_profile.set_title("Detailed Uniformity Profiles")
        ax_profile.set_xlabel("Angle (deg)")
        ax_profile.set_ylabel("Pixel Value (HU)")
        ax_profile.grid(True, alpha=0.3)
        ax_profile.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=9)

        ax_smooth.set_title("Smoothed Profiles (Savitzky-Golay)")
        ax_smooth.set_xlabel("Angle (deg)")
        ax_smooth.set_ylabel("Pixel Value (HU)")
        ax_smooth.grid(True, alpha=0.3)
        ax_smooth.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)

        ax_resid.set_title("Residuals vs Angle (Value - Mean)")
        ax_resid.set_xlabel("Angle (deg)")
        ax_resid.set_ylabel("Residual (HU)")
        ax_resid.grid(True, alpha=0.3)
        ax_resid.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)

        ax_hist.set_title("Profile Value Histograms")
        ax_hist.set_xlabel("Pixel Value (HU)")
        ax_hist.set_ylabel("Count")
        ax_hist.grid(True, alpha=0.3)
        ax_hist.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)

        all_profiles = self._get_all_profiles()
        mean_points_all = [float(np.mean(p["values"])) for p in all_profiles]
        std_points_all = [float(np.std(p["values"])) for p in all_profiles]
        radius_points_all = [float(p["radius_mm"]) for p in all_profiles]

        ax_mean.errorbar(radius_points_all, mean_points_all, yerr=std_points_all, fmt='o-', color='black', ecolor='gray', capsize=4, label='Mean (HU)')
        ax_mean.set_title("Mean and Std vs Radius")
        ax_mean.set_xlabel("Radius (mm)")
        ax_mean.set_ylabel("Mean Pixel Value (HU)")
        ax_mean.grid(True, alpha=0.3)
        ax_std = ax_mean.twinx()
        ax_std.plot(radius_points_all, std_points_all, 's--', color='tab:blue', label='Std Dev (HU)')
        ax_std.set_ylabel("Std Dev (HU)")
        lines_mean, labels_mean = ax_mean.get_legend_handles_labels()
        lines_std, labels_std = ax_std.get_legend_handles_labels()
        ax_mean.legend(lines_mean + lines_std, labels_mean + labels_std, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=9)

        image = getattr(self.analyzer, "image", None)
        center = getattr(self.analyzer, "center", None)
        if image is None or center is None:
            raise ValueError("Analyzer image and center are required for overlay plot")

        ax_img.imshow(image, cmap="gray")
        ax_img.set_title("Sampling Rings")
        ax_img.axis("off")

        cx, cy = center
        for idx, profile in enumerate(profiles):
            radius_mm = profile["radius_mm"]
            radius_px = radius_mm / float(getattr(self.analyzer, "pixel_spacing", 1.0))
            theta = np.linspace(0, 2 * np.pi, 200)
            circle_x = cx + radius_px * np.cos(theta)
            circle_y = cy + radius_px * np.sin(theta)
            ax_img.plot(circle_x, circle_y, color=colors[idx % len(colors)], linewidth=1.2)

        ax_img.plot(cx, cy, "r+", markersize=10, markeredgewidth=2)
        fig.tight_layout()
        return fig
