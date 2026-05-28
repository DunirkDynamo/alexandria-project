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

    def _angles_from_profile(self, profile: Dict[str, Any]) -> np.ndarray:
        if "angles_deg" in profile:
            return np.array(profile["angles_deg"], dtype=float)

        values = np.array(profile.get("values", []), dtype=float)
        start = float(profile.get("start_angle_deg", 0.0))
        increment = float(profile.get("angle_increment_deg", 0.0))
        if values.size == 0:
            return np.array([], dtype=float)
        return start + increment * np.arange(values.size, dtype=float)

    def _get_profiles(self) -> List[Dict[str, Any]]:
        plot_radii = set(round(float(r), 3) for r in getattr(self.analyzer, "radii_mm", self.analyzer.results.get("plot_radii_mm", [])))
        if getattr(self.analyzer, "profile_data", None):
            return [profile for profile in self.analyzer.profile_data if round(float(profile.get("radius_mm", 0.0)), 3) in plot_radii]
        profiles = []
        for profile in self.analyzer.results.get("profiles", []):
            if round(float(profile.get("radius_mm", 0.0)), 3) not in plot_radii:
                continue
            profiles.append({
                "radius_mm": profile["radius_mm"],
                "angles_deg": self._angles_from_profile(profile),
                "values": np.array(profile["values"], dtype=float),
            })
        return profiles

    def _get_all_profiles(self) -> List[Dict[str, Any]]:
        if getattr(self.analyzer, "profile_data", None):
            return self.analyzer.profile_data
        profiles = []
        for profile in self.analyzer.results.get("profiles", []):
            profiles.append({
                "radius_mm": profile["radius_mm"],
                "angles_deg": self._angles_from_profile(profile),
                "values": np.array(profile["values"], dtype=float),
            })
        return profiles

    def plot(self, bins: int = 25, figsize: tuple = (14, 20)) -> plt.Figure:
        profiles = self._get_profiles()
        if not profiles:
            raise ValueError("No profile data available for plotting")

        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(8, 3, height_ratios=[1.1, 1, 1, 1, 1, 1, 1.1, 0.25], width_ratios=[1, 1, 1])
        ax_rings = fig.add_subplot(gs[0, :2])
        ax_profile = fig.add_subplot(gs[1, :2])
        ax_smooth = fig.add_subplot(gs[2, :2])
        ax_resid = fig.add_subplot(gs[3, :2])
        ax_hist = fig.add_subplot(gs[4, 0])
        ax_fit = fig.add_subplot(gs[4, 1], sharey=ax_hist)
        ax_mean = fig.add_subplot(gs[5, :2])
        ax_img = fig.add_subplot(gs[6, 0])
        ax_cross = fig.add_subplot(gs[6, 1])
        ax_cross_legend = fig.add_subplot(gs[7, :2])
        colors = plt.cm.tab10(np.linspace(0, 1, 7))

        mean_points = []
        std_points = []
        radius_points = []

        fit_lines = []
        fit_labels = []
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
            # Histogram
            ax_hist.hist(values, bins=bins, histtype='step', color=color, linewidth=2, label=f"r={radius:.1f}mm")
            # Gaussian fit
            std_val = float(np.std(values))
            if std_val > 0.0:
                x_min = float(np.min(values))
                x_max = float(np.max(values))
                if x_max > x_min:
                    x_vals = np.linspace(x_min, x_max, 200)
                    bin_width = (x_max - x_min) / float(bins)
                    pdf = (1.0 / (std_val * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x_vals - float(np.mean(values))) / std_val) ** 2)
                    line, = ax_fit.plot(x_vals, pdf * len(values) * bin_width, color=color, linestyle='--', linewidth=1.5, label=f"r={radius:.1f}mm")
                    fit_lines.append(line)
                    fit_labels.append(f"r={radius:.1f}mm")
            radius_points.append(radius)
            mean_points.append(mean_val)
            std_points.append(std_val)

        ax_profile.set_title("Detailed Uniformity Profiles")
        ax_profile.set_xlabel("Angle (deg)")
        ax_profile.set_ylabel("Pixel Value (HU)")
        ax_profile.grid(True, alpha=0.3)
        ax_profile.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=9)

        # Show Savitzky-Golay parameters in the plot title
        sg_window = min(17, len(profiles[0]["values"]) - 1 if len(profiles[0]["values"]) % 2 == 0 else len(profiles[0]["values"])) if profiles else 17
        if sg_window % 2 == 0:
            sg_window -= 1
        sg_poly = 3
        ax_smooth.set_title(f"Smoothed Profiles (Savitzky-Golay, window={sg_window}, polyorder={sg_poly})")
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
        # No legend on left

        ax_fit.set_title("Gaussian Fits")
        ax_fit.set_xlabel("Pixel Value (HU)")
        ax_fit.grid(True, alpha=0.3)
        # Legend to the right of the fit plot
        if fit_lines:
            ax_fit.legend(fit_lines, fit_labels, loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)


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

        # Top: Sampling rings image
        ax_rings.imshow(image, cmap="gray")
        ax_rings.set_title("Sampling Rings")
        ax_rings.axis("off")
        cx, cy = center
        for idx, profile in enumerate(profiles):
            radius_mm = profile["radius_mm"]
            radius_px = radius_mm / float(getattr(self.analyzer, "pixel_spacing", 1.0))
            theta = np.linspace(0, 2 * np.pi, 200)
            circle_x = cx + radius_px * np.cos(theta)
            circle_y = cy + radius_px * np.sin(theta)
            ax_rings.plot(circle_x, circle_y, color=colors[idx % len(colors)], linewidth=1.2)
        ax_rings.plot(cx, cy, "r+", markersize=10, markeredgewidth=2)

        # --- Cross-section profiles at 0, 30, ..., 150 deg ---
        cross_angles_deg_math = np.arange(0, 180, 30)
        cross_angles_deg_img = -cross_angles_deg_math  # Convert to image convention (clockwise positive)
        n_cross = len(cross_angles_deg_math)
        cross_colors = plt.cm.tab10(np.linspace(0, 1, n_cross))
        profile_len = 400
        half_len = profile_len // 2
        pixel_spacing = float(getattr(self.analyzer, "pixel_spacing", 1.0))

        cross_profiles = []
        dists = np.arange(-half_len, half_len, 1)
        for i, angle_deg_img in enumerate(cross_angles_deg_img):
            theta = np.deg2rad(angle_deg_img)
            dx = np.cos(theta)
            dy = np.sin(theta)
            xs = cx + dists * dx
            ys = cy + dists * dy
            xs_clip = np.clip(xs, 0, image.shape[1] - 1)
            ys_clip = np.clip(ys, 0, image.shape[0] - 1)
            profile = image[ys_clip.astype(int), xs_clip.astype(int)]
            cross_profiles.append(profile)

        # Plot cross-section profiles (right column)
        ax_cross.set_title("Center Cross-Sections (0° to 150°)")
        cross_lines = []
        cross_labels = []
        for i, (angle_deg, profile) in enumerate(zip(cross_angles_deg_math, cross_profiles)):
            line, = ax_cross.plot(dists, profile, color=cross_colors[i], label=f"{int(angle_deg)}°")
            cross_lines.append(line)
            cross_labels.append(f"{int(angle_deg)}°")
        ax_cross.set_xlabel("Distance from Center (px)")
        ax_cross.set_ylabel("Pixel Value (HU)")
        ax_cross.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
        ax_cross.grid(True, alpha=0.3)

        # Place legend between image and plot
        ax_cross_legend.axis('off')
        ax_cross_legend.legend(cross_lines, cross_labels, loc='center', ncol=4, fontsize=10, frameon=False)

        # Show cross-section lines on image (left column)
        ax_img.imshow(image, cmap="gray")
        ax_img.set_title("Cross-Sections Overlay")
        ax_img.axis("off")
        for i, angle_deg_img in enumerate(cross_angles_deg_img):
            theta = np.deg2rad(angle_deg_img)
            x0 = cx - half_len * np.cos(theta)
            y0 = cy - half_len * np.sin(theta)
            x1 = cx + half_len * np.cos(theta)
            y1 = cy + half_len * np.sin(theta)
            ax_img.plot([x0, x1], [y0, y1], color=cross_colors[i], linewidth=2, label=f"{int(cross_angles_deg_math[i])}°")
        ax_img.plot(cx, cy, "r+", markersize=8, markeredgewidth=1.5)

        fig.tight_layout()
        return fig

        fig.tight_layout()
        return fig
