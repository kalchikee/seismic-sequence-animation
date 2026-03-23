"""
animation.py
------------
Creates an MP4 animation of the 2019 Ridgecrest seismic sequence.

Each frame shows cumulative seismicity up to that date with:
  - Circle size proportional to magnitude (area ∝ 10^M)
  - Color indicating depth (shallow = red, deep = blue)
  - Date stamp and running event count
  - Fault traces and region basemap

Requirements: ffmpeg must be installed (conda install ffmpeg or brew install ffmpeg)
"""

import logging
from pathlib import Path

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

matplotlib.use("Agg")
log = logging.getLogger(__name__)

# Style constants
DEPTH_CMAP = "RdYlBu_r"
DEPTH_VMIN, DEPTH_VMAX = 0.0, 25.0
FPS = 10
DPI = 150
FIGURE_SIZE = (12, 10)

STUDY_BOUNDS = {
    "xlim": (-118.8, -116.2),
    "ylim": (34.8, 36.7),
}


def magnitude_to_size(magnitude: np.ndarray, base: float = 4.0) -> np.ndarray:
    """Scale circle area by 10^M for visual magnitude representation."""
    return base * 10 ** (0.8 * (magnitude - 1.5))


def create_animation(
    catalog: gpd.GeoDataFrame,
    faults: gpd.GeoDataFrame,
    output_path: Path = Path("results/ridgecrest_sequence.mp4"),
    time_col: str = "days_after_mainshock",
    frames_per_day: int = 1,
) -> None:
    """
    Build and save the seismicity animation.

    Parameters
    ----------
    catalog : GeoDataFrame
        Earthquake catalog with time, magnitude, depth_km, latitude, longitude.
    faults : GeoDataFrame
        Fault traces to overlay.
    output_path : Path
        Output MP4 file path.
    time_col : str
        Column name for days after mainshock.
    frames_per_day : int
        Number of animation frames per calendar day.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Sort catalog chronologically
    catalog = catalog.sort_values(time_col).reset_index(drop=True)

    t_min = catalog[time_col].min()
    t_max = catalog[time_col].max()
    frame_times = np.arange(t_min, t_max + 1.0 / frames_per_day, 1.0 / frames_per_day)
    n_frames = len(frame_times)

    log.info(f"Creating animation: {n_frames} frames, {FPS} fps → {output_path}")

    # Pre-compute depth colormap
    norm = Normalize(vmin=DEPTH_VMIN, vmax=DEPTH_VMAX)
    cmap = plt.cm.get_cmap(DEPTH_CMAP)
    sm = ScalarMappable(norm=norm, cmap=cmap)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI)
    fig.patch.set_facecolor("#1a1a2e")
    ax.set_facecolor("#1a1a2e")

    # ── Static elements ────────────────────────────────────────────────────────

    # Fault traces
    if faults is not None and len(faults) > 0:
        faults_reproj = faults.to_crs("EPSG:4326") if faults.crs else faults
        faults_reproj.plot(
            ax=ax, color="#FFD700", linewidth=0.6, alpha=0.7, zorder=2
        )

    # Mainshock star
    ax.plot(
        -117.599, 35.770, "*", color="white", markersize=18,
        markeredgecolor="black", markeredgewidth=1.0, zorder=10,
        label="M7.1 Mainshock"
    )

    ax.set_xlim(STUDY_BOUNDS["xlim"])
    ax.set_ylim(STUDY_BOUNDS["ylim"])
    ax.set_xlabel("Longitude", color="white", fontsize=11)
    ax.set_ylabel("Latitude", color="white", fontsize=11)
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#444")

    # Colorbar
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Depth (km)", color="white", fontsize=10)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    # Magnitude legend
    for mag, label in [(2.0, "M2"), (3.0, "M3"), (4.0, "M4"), (5.0, "M5")]:
        ax.scatter(
            [], [], s=magnitude_to_size(np.array([mag])),
            c="gray", alpha=0.7, label=label, edgecolors="none"
        )

    legend = ax.legend(
        loc="lower left", fontsize=9, framealpha=0.3,
        facecolor="#1a1a2e", labelcolor="white"
    )

    # Title and date stamp
    title = ax.set_title(
        "2019 Ridgecrest Earthquake Sequence", color="white",
        fontsize=14, fontweight="bold", pad=10
    )
    date_text = ax.text(
        0.02, 0.96, "", transform=ax.transAxes,
        color="white", fontsize=11, va="top",
        fontfamily="monospace"
    )
    count_text = ax.text(
        0.02, 0.91, "", transform=ax.transAxes,
        color="#aaaaaa", fontsize=9, va="top"
    )

    # Dynamic scatter (starts empty)
    sc = ax.scatter(
        [], [], c=[], cmap=DEPTH_CMAP, norm=norm,
        s=[], alpha=0.65, edgecolors="none", zorder=5
    )

    # ── Animation update function ──────────────────────────────────────────────

    def update(frame_idx: int):
        t = frame_times[frame_idx]
        visible = catalog[catalog[time_col] <= t]

        if len(visible) > 0:
            lons = visible["longitude"].values
            lats = visible["latitude"].values
            depths = np.clip(visible["depth_km"].values, DEPTH_VMIN, DEPTH_VMAX)
            sizes = magnitude_to_size(visible["magnitude"].values)

            sc.set_offsets(np.column_stack([lons, lats]))
            sc.set_array(depths)
            sc.set_sizes(sizes)
        else:
            sc.set_offsets(np.empty((0, 2)))

        # Date label
        date_label = (
            pd.Timestamp("2019-07-06") + pd.Timedelta(days=t)
        ).strftime("%Y-%m-%d")
        date_text.set_text(f"Date: {date_label}")
        count_text.set_text(f"Cumulative events (M≥1.5): {len(visible):,}")

        return sc, date_text, count_text

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames, interval=1000 / FPS,
        blit=True, repeat=False
    )

    writer = animation.FFMpegWriter(fps=FPS, bitrate=2000)
    ani.save(str(output_path), writer=writer)
    plt.close(fig)
    log.info(f"Animation saved → {output_path}")


def export_key_frames(
    catalog: gpd.GeoDataFrame,
    faults: gpd.GeoDataFrame,
    output_dir: Path = Path("results/frames"),
    days: list[float] | None = None,
) -> None:
    """Export static PNG frames at key time steps for portfolio presentation."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if days is None:
        days = [1.0, 7.0, 30.0, 90.0, 180.0]

    norm = Normalize(vmin=DEPTH_VMIN, vmax=DEPTH_VMAX)

    for day in days:
        visible = catalog[catalog["days_after_mainshock"] <= day]
        fig, ax = plt.subplots(figsize=(10, 8), dpi=DPI)
        fig.patch.set_facecolor("#1a1a2e")
        ax.set_facecolor("#1a1a2e")

        if faults is not None and len(faults) > 0:
            faults.to_crs("EPSG:4326").plot(ax=ax, color="#FFD700", linewidth=0.7, alpha=0.7)

        if len(visible) > 0:
            sc = ax.scatter(
                visible["longitude"], visible["latitude"],
                c=visible["depth_km"].clip(DEPTH_VMIN, DEPTH_VMAX),
                s=magnitude_to_size(visible["magnitude"].values),
                cmap=DEPTH_CMAP, norm=norm, alpha=0.65,
                edgecolors="none", zorder=5
            )
            fig.colorbar(sc, ax=ax, label="Depth (km)", shrink=0.6)

        ax.plot(-117.599, 35.770, "*", color="white", markersize=16, zorder=10)
        ax.set_xlim(STUDY_BOUNDS["xlim"])
        ax.set_ylim(STUDY_BOUNDS["ylim"])
        ax.set_title(
            f"Ridgecrest Sequence — Day {day:.0f} ({len(visible):,} events)",
            color="white", fontsize=13
        )
        ax.tick_params(colors="white")

        fname = output_dir / f"ridgecrest_day_{int(day):03d}.png"
        fig.savefig(fname, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)
        log.info(f"Frame saved → {fname}")
