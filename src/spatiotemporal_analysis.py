"""
spatiotemporal_analysis.py
--------------------------
Spatiotemporal analysis of the 2019 Ridgecrest aftershock sequence:
  - Omori law fitting (K, c, p parameters with bootstrap uncertainties)
  - Along-strike / along-dip distance calculations from mainshock fault plane
  - Seismicity migration analysis (afterslip front vs. fluid diffusion)
  - B-value spatial mapping using Aki maximum likelihood estimator

References:
  Utsu (1961) Omori law; Aki (1965) b-value MLE
"""

import logging
from pathlib import Path
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm

log = logging.getLogger(__name__)

# Mainshock fault plane parameters (Goldberg et al. 2020)
MAINSHOCK_LAT = 35.770
MAINSHOCK_LON = -117.599
FAULT_STRIKE_DEG = 322.0  # clockwise from north
FAULT_DIP_DEG = 81.0      # from horizontal


# ── Coordinate helpers ─────────────────────────────────────────────────────────

def latlon_to_utm(lat: np.ndarray, lon: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Convert geographic coordinates to UTM Zone 11N (meters)."""
    import pyproj
    proj = pyproj.Proj(proj="utm", zone=11, ellps="WGS84")
    x, y = proj(lon, lat)
    return x, y


def compute_fault_distances(
    gdf: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Project each hypocenter onto the mainshock fault plane coordinate system.

    Adds columns:
      along_strike_km : distance along fault strike (positive = NW)
      along_dip_km    : distance down-dip (positive = deeper)
      fault_normal_km : distance perpendicular to fault plane

    Parameters
    ----------
    gdf : GeoDataFrame with latitude, longitude, depth_km columns

    Returns
    -------
    GeoDataFrame with added distance columns
    """
    gdf = gdf.copy()

    # UTM coordinates
    x, y = latlon_to_utm(gdf["latitude"].values, gdf["longitude"].values)
    x0, y0 = latlon_to_utm(
        np.array([MAINSHOCK_LAT]), np.array([MAINSHOCK_LON])
    )

    dx = x - x0[0]
    dy = y - y0[0]
    dz = (gdf["depth_km"].values - 8.0) * 1000.0  # depth relative to mainshock (m)

    # Strike and dip unit vectors
    strike_rad = np.radians(FAULT_STRIKE_DEG)
    dip_rad = np.radians(FAULT_DIP_DEG)

    # Along-strike unit vector (horizontal, in direction of strike)
    s = np.array([np.sin(strike_rad), np.cos(strike_rad), 0.0])

    # Down-dip unit vector
    d = np.array([
        np.cos(strike_rad) * np.cos(dip_rad),
        -np.sin(strike_rad) * np.cos(dip_rad),
        -np.sin(dip_rad),
    ])

    # Fault-normal unit vector
    n = np.cross(s, d)

    coords = np.column_stack([dx, dy, dz])

    gdf["along_strike_km"] = coords @ s / 1000.0
    gdf["along_dip_km"] = coords @ d / 1000.0
    gdf["fault_normal_km"] = coords @ n / 1000.0

    log.info(
        f"Strike range: {gdf['along_strike_km'].min():.1f} to "
        f"{gdf['along_strike_km'].max():.1f} km"
    )
    return gdf


# ── Omori Law ──────────────────────────────────────────────────────────────────

def omori_rate(t: np.ndarray, K: float, c: float, p: float) -> np.ndarray:
    """Modified Omori aftershock rate: n(t) = K / (t + c)^p"""
    return K / (t + c) ** p


def fit_omori_law(
    gdf: gpd.GeoDataFrame,
    time_col: str = "days_after_mainshock",
    dt_bin: float = 1.0,
    n_bootstrap: int = 500,
) -> dict:
    """
    Fit modified Omori law to aftershock rate time series.

    Returns
    -------
    dict with keys: K, c, p, K_std, c_std, p_std, t_bins, rate_observed, rate_fitted
    """
    # Keep only aftershocks (t > 0)
    after = gdf[gdf[time_col] > 0].copy()
    t_max = after[time_col].max()

    bins = np.arange(0, t_max + dt_bin, dt_bin)
    counts, _ = np.histogram(after[time_col], bins=bins)
    t_centers = 0.5 * (bins[:-1] + bins[1:])

    # Mask zero-count bins for fitting
    mask = counts > 0
    t_fit = t_centers[mask]
    n_fit = counts[mask].astype(float)

    p0 = [100.0, 0.1, 1.0]
    bounds = ([0, 1e-5, 0.5], [1e6, 10.0, 2.0])

    try:
        popt, pcov = curve_fit(
            omori_rate, t_fit, n_fit, p0=p0, bounds=bounds, maxfev=10000
        )
        K, c, p = popt
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError:
        log.warning("Omori curve_fit failed; using initial estimates")
        K, c, p = p0
        perr = np.array([np.nan, np.nan, np.nan])

    # Bootstrap uncertainty
    boot_params = []
    rng = np.random.default_rng(42)
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(t_fit), size=len(t_fit))
        try:
            popt_b, _ = curve_fit(
                omori_rate, t_fit[idx], n_fit[idx], p0=[K, c, p],
                bounds=bounds, maxfev=5000
            )
            boot_params.append(popt_b)
        except RuntimeError:
            pass

    if boot_params:
        boot_arr = np.array(boot_params)
        perr = boot_arr.std(axis=0)

    rate_fitted = omori_rate(t_centers, K, c, p)

    log.info(f"Omori fit: K={K:.1f}±{perr[0]:.1f}, c={c:.4f}±{perr[1]:.4f}, p={p:.3f}±{perr[2]:.3f}")

    return {
        "K": K, "c": c, "p": p,
        "K_std": perr[0], "c_std": perr[1], "p_std": perr[2],
        "t_bins": t_centers,
        "rate_observed": counts,
        "rate_fitted": rate_fitted,
    }


# ── B-Value Analysis ───────────────────────────────────────────────────────────

def aki_bvalue(magnitudes: np.ndarray, mc: float) -> Tuple[float, float]:
    """
    Aki (1965) maximum likelihood b-value estimator.

    b = log10(e) / (M_mean - Mc)

    Returns
    -------
    b : float
        Gutenberg-Richter b-value
    sigma_b : float
        Standard error (Shi & Bolt 1982)
    """
    mags = magnitudes[magnitudes >= mc]
    if len(mags) < 10:
        return np.nan, np.nan
    m_mean = mags.mean()
    b = np.log10(np.e) / (m_mean - mc + 0.05)  # 0.05 = bin width/2
    sigma_b = b / np.sqrt(len(mags))
    return b, sigma_b


def completeness_magnitude(magnitudes: np.ndarray, bin_width: float = 0.1) -> float:
    """
    Estimate completeness magnitude (Mc) using maximum curvature method.

    Mc is the magnitude bin with the highest frequency of occurrence.
    """
    bins = np.arange(magnitudes.min(), magnitudes.max() + bin_width, bin_width)
    counts, edges = np.histogram(magnitudes, bins=bins)
    mc = edges[np.argmax(counts)] + bin_width / 2
    return round(mc, 1)


def spatial_bvalue_map(
    gdf: gpd.GeoDataFrame,
    bin_deg: float = 0.1,
    min_events: int = 30,
) -> pd.DataFrame:
    """
    Compute b-values in spatial bins.

    Returns DataFrame with columns: lat_center, lon_center, b, sigma_b, n_events, mc
    """
    lats = np.arange(
        gdf["latitude"].min(), gdf["latitude"].max() + bin_deg, bin_deg
    )
    lons = np.arange(
        gdf["longitude"].min(), gdf["longitude"].max() + bin_deg, bin_deg
    )

    records = []
    for lat in lats[:-1]:
        for lon in lons[:-1]:
            subset = gdf[
                (gdf["latitude"] >= lat) & (gdf["latitude"] < lat + bin_deg) &
                (gdf["longitude"] >= lon) & (gdf["longitude"] < lon + bin_deg)
            ]
            if len(subset) < min_events:
                continue
            mags = subset["magnitude"].values
            mc = completeness_magnitude(mags)
            b, sigma_b = aki_bvalue(mags, mc)
            records.append({
                "lat_center": lat + bin_deg / 2,
                "lon_center": lon + bin_deg / 2,
                "b": b,
                "sigma_b": sigma_b,
                "n_events": len(subset),
                "mc": mc,
            })

    df = pd.DataFrame(records)
    log.info(f"Computed b-values in {len(df)} spatial bins")
    return df


# ── Migration Analysis ─────────────────────────────────────────────────────────

def migration_analysis(gdf: gpd.GeoDataFrame) -> dict:
    """
    Test whether seismicity migrates along strike over time.

    Fits:
      1. Linear migration: r = v * t
      2. Diffusive migration: r = sqrt(4 * D * t)

    Returns
    -------
    dict with migration velocity, diffusivity, and R² values
    """
    after = gdf[
        (gdf["days_after_mainshock"] > 0) &
        (gdf["days_after_mainshock"] < 30)
    ].copy()

    t = after["days_after_mainshock"].values
    r = np.abs(after["along_strike_km"].values)

    # Linear fit
    A_lin = np.column_stack([t, np.ones_like(t)])
    result_lin = np.linalg.lstsq(A_lin, r, rcond=None)
    v_km_per_day = result_lin[0][0]
    r_lin = A_lin @ result_lin[0]
    ss_res = np.sum((r - r_lin) ** 2)
    ss_tot = np.sum((r - r.mean()) ** 2)
    r2_lin = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # Diffusive fit: r = sqrt(4Dt)
    sqrt_t = np.sqrt(t)
    D_fit = np.linalg.lstsq(
        sqrt_t.reshape(-1, 1), r, rcond=None
    )[0][0]
    D_km2_per_day = (D_fit ** 2) / 4.0
    r_diff = D_fit * sqrt_t
    ss_res_d = np.sum((r - r_diff) ** 2)
    r2_diff = 1 - ss_res_d / ss_tot if ss_tot > 0 else 0

    preferred = "linear" if r2_lin > r2_diff else "diffusive"
    log.info(
        f"Migration: v={v_km_per_day:.2f} km/day (R²={r2_lin:.3f}); "
        f"D={D_km2_per_day:.2f} km²/day (R²={r2_diff:.3f}); "
        f"preferred={preferred}"
    )

    return {
        "velocity_km_per_day": v_km_per_day,
        "r2_linear": r2_lin,
        "diffusivity_km2_per_day": D_km2_per_day,
        "r2_diffusive": r2_diff,
        "preferred_model": preferred,
    }
