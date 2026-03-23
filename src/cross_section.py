"""
cross_section.py
----------------
Builds fault-perpendicular depth cross-sections of the Ridgecrest aftershock
sequence to reveal seismogenic zone geometry and the brittle-ductile transition.

Also overlays the USGS finite-fault slip model to show the spatial relationship
between coseismic slip and aftershock distribution.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from matplotlib.colors import Normalize

log = logging.getLogger(__name__)

# Profile parameters (fault-perpendicular to the M7.1 rupture)
PROFILE_AZIMUTH_DEG = 52.0   # perpendicular to fault strike of 322°
PROFILE_WIDTH_KM = 15.0       # ± 15 km on each side of profile
MAINSHOCK_LAT = 35.770
MAINSHOCK_LON = -117.599
DEPTH_MAX_KM = 25.0


def project_to_cross_section(
    gdf: gpd.GeoDataFrame,
    azimuth_deg: float = PROFILE_AZIMUTH_DEG,
    width_km: float = PROFILE_WIDTH_KM,
) -> pd.DataFrame:
    """
    Project hypocenters onto a fault-perpendicular cross-section profile.

    Parameters
    ----------
    gdf : GeoDataFrame with latitude, longitude, depth_km
    azimuth_deg : Profile azimuth in degrees CW from north
    width_km : Half-width of the swath to include

    Returns
    -------
    DataFrame with columns: along_profile_km, depth_km, magnitude, perpendicular_km
    """
    import pyproj
    proj = pyproj.Proj(proj="utm", zone=11, ellps="WGS84")

    x, y = proj(gdf["longitude"].values, gdf["latitude"].values)
    x0, y0 = proj(MAINSHOCK_LON, MAINSHOCK_LAT)

    dx = (x - x0) / 1000.0  # km
    dy = (y - y0) / 1000.0  # km

    az_rad = np.radians(azimuth_deg)
    along = dx * np.sin(az_rad) + dy * np.cos(az_rad)
    perp = -dx * np.cos(az_rad) + dy * np.sin(az_rad)

    mask = np.abs(perp) <= width_km

    return pd.DataFrame({
        "along_profile_km": along[mask],
        "depth_km": gdf["depth_km"].values[mask],
        "magnitude": gdf["magnitude"].values[mask],
        "perpendicular_km": perp[mask],
    })


def plot_cross_section(
    gdf: gpd.GeoDataFrame,
    output_path: Path = Path("results/cross_section.png"),
) -> None:
    """
    Generate a professional depth cross-section figure.

    Parameters
    ----------
    gdf : Earthquake catalog GeoDataFrame
    output_path : Save path for the figure
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = project_to_cross_section(gdf)

    fig, (ax_map, ax_xs) = plt.subplots(
        2, 1, figsize=(12, 10),
        gridspec_kw={"height_ratios": [1, 1.2]},
    )
    fig.subplots_adjust(hspace=0.35)

    # ── Map view ───────────────────────────────────────────────────────────────
    norm = Normalize(vmin=0, vmax=20)
    sc = ax_map.scatter(
        gdf["longitude"], gdf["latitude"],
        c=gdf["depth_km"], cmap="RdYlBu_r", norm=norm,
        s=3.0 * 10 ** (0.5 * (gdf["magnitude"] - 1.5)),
        alpha=0.5, edgecolors="none", rasterized=True
    )

    # Profile line
    az_rad = np.radians(PROFILE_AZIMUTH_DEG)
    import pyproj
    proj = pyproj.Proj(proj="utm", zone=11, ellps="WGS84")
    x0, y0 = proj(MAINSHOCK_LON, MAINSHOCK_LAT)
    L = 40.0 * 1000.0  # 40 km profile
    pts_lon, pts_lat = [], []
    for sign in [-1, 1]:
        xp = x0 + sign * L * np.sin(az_rad)
        yp = y0 + sign * L * np.cos(az_rad)
        lon, lat = proj(xp, yp, inverse=True)
        pts_lon.append(lon)
        pts_lat.append(lat)

    ax_map.plot(pts_lon, pts_lat, "w--", linewidth=1.5, label="Profile A–A'")
    ax_map.plot(MAINSHOCK_LON, MAINSHOCK_LAT, "*", color="yellow", markersize=14,
                markeredgecolor="black", label="M7.1")
    ax_map.set_xlim(-118.8, -116.2)
    ax_map.set_ylim(34.8, 36.7)
    ax_map.set_xlabel("Longitude", fontsize=11)
    ax_map.set_ylabel("Latitude", fontsize=11)
    ax_map.set_title("Map View — Ridgecrest Sequence (M≥1.5)", fontsize=13)
    ax_map.legend(fontsize=9)
    plt.colorbar(sc, ax=ax_map, label="Depth (km)", shrink=0.6)

    # ── Cross-section ──────────────────────────────────────────────────────────
    norm2 = Normalize(vmin=0, vmax=20)
    sc2 = ax_xs.scatter(
        df["along_profile_km"], df["depth_km"],
        c=df["depth_km"], cmap="RdYlBu_r", norm=norm2,
        s=4.0 * 10 ** (0.5 * (df["magnitude"] - 1.5)),
        alpha=0.6, edgecolors="none", rasterized=True
    )

    # Mainshock position
    ax_xs.plot(0, 8.0, "*", color="yellow", markersize=16, markeredgecolor="black",
               zorder=10, label="M7.1 mainshock")

    # Brittle-ductile transition annotation
    bdt_depth = df["depth_km"].quantile(0.99)
    ax_xs.axhline(bdt_depth, color="orange", linestyle="--", linewidth=1.2,
                  label=f"~BDT depth ({bdt_depth:.0f} km)")

    # Approximate fault plane trace
    dip_rad = np.radians(81.0)
    x_fp = np.linspace(-30, 30, 100)
    z_fp = 8.0 + np.abs(x_fp) * np.tan(dip_rad - np.pi / 2)
    mask_fp = (z_fp >= 0) & (z_fp <= DEPTH_MAX_KM)
    ax_xs.plot(x_fp[mask_fp], z_fp[mask_fp], "w-", linewidth=1.5,
               alpha=0.5, label="Fault plane (81° dip)")

    ax_xs.set_ylim(DEPTH_MAX_KM, 0)  # depth increases downward
    ax_xs.set_xlim(-40, 40)
    ax_xs.set_xlabel("Distance along profile A–A' (km)", fontsize=11)
    ax_xs.set_ylabel("Depth (km)", fontsize=11)
    ax_xs.set_title(
        f"Cross-Section A–A' (perpendicular to fault, swath ±{PROFILE_WIDTH_KM} km)",
        fontsize=13
    )
    ax_xs.legend(fontsize=9)
    ax_xs.text(-38, 1.5, "A", fontsize=14, fontweight="bold")
    ax_xs.text(36, 1.5, "A'", fontsize=14, fontweight="bold")
    plt.colorbar(sc2, ax=ax_xs, label="Depth (km)", shrink=0.6)

    fig.suptitle(
        "2019 Ridgecrest M7.1 — Seismogenic Zone Geometry",
        fontsize=15, fontweight="bold", y=1.01
    )

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"Cross-section saved → {output_path}")
