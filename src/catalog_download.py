"""
catalog_download.py
-------------------
Downloads the 2019 Ridgecrest earthquake sequence from the USGS FDSN web service
and stores the catalog as a GeoPackage for use in QGIS and Python analysis.

Data source: USGS Earthquake Hazards Program
API: https://earthquake.usgs.gov/fdsnws/event/1/query
"""

import logging
import time
from pathlib import Path
from typing import Optional

import geopandas as gpd
import pandas as pd
import requests
from shapely.geometry import Point

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)

# ── Configuration ──────────────────────────────────────────────────────────────
USGS_FDSN_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
QFFDB_URL = (
    "https://earthquake.usgs.gov/static/lfs/nshm/qfaults/Qfaults_GIS_2021_new.zip"
)

STUDY_REGION = {
    "minlatitude": 35.0,
    "maxlatitude": 36.5,
    "minlongitude": -118.5,
    "maxlongitude": -116.5,
}

MAINSHOCK = {
    "event_id": "ci38457511",
    "magnitude": 7.1,
    "latitude": 35.770,
    "longitude": -117.599,
    "depth_km": 8.0,
    "time": "2019-07-06T03:19:53",
    "strike": 322,
    "dip": 81,
}


def download_catalog(
    starttime: str = "2019-07-01",
    endtime: str = "2020-01-01",
    minmagnitude: float = 1.5,
    output_dir: Path = Path("data"),
    max_retries: int = 3,
) -> gpd.GeoDataFrame:
    """
    Query the USGS FDSN earthquake catalog and return a GeoDataFrame.

    Parameters
    ----------
    starttime : str
        ISO 8601 start date.
    endtime : str
        ISO 8601 end date.
    minmagnitude : float
        Minimum magnitude threshold.
    output_dir : Path
        Directory to save GeoPackage and CSV.
    max_retries : int
        Number of HTTP retry attempts.

    Returns
    -------
    GeoDataFrame
        Earthquake catalog with geometry (Point) and attributes.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "format": "geojson",
        "starttime": starttime,
        "endtime": endtime,
        "minmagnitude": minmagnitude,
        **STUDY_REGION,
        "orderby": "time",
        "limit": 20000,
    }

    log.info(
        f"Querying USGS FDSN: M≥{minmagnitude} events, "
        f"{starttime} to {endtime}, region {STUDY_REGION}"
    )

    response = None
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(USGS_FDSN_URL, params=params, timeout=60)
            response.raise_for_status()
            break
        except requests.RequestException as exc:
            log.warning(f"Attempt {attempt}/{max_retries} failed: {exc}")
            if attempt < max_retries:
                time.sleep(5 * attempt)
            else:
                raise

    data = response.json()
    features = data["features"]
    log.info(f"Retrieved {len(features):,} events from USGS ComCat")

    records = []
    for feat in features:
        props = feat["properties"]
        coords = feat["geometry"]["coordinates"]
        records.append(
            {
                "event_id": feat["id"],
                "time": pd.to_datetime(props["time"], unit="ms", utc=True),
                "magnitude": props["mag"],
                "mag_type": props.get("magType", ""),
                "depth_km": coords[2],
                "longitude": coords[0],
                "latitude": coords[1],
                "place": props.get("place", ""),
                "status": props.get("status", ""),
                "net": props.get("net", ""),
            }
        )

    df = pd.DataFrame(records)
    df["days_after_mainshock"] = (
        df["time"] - pd.Timestamp("2019-07-06 03:19:53", tz="UTC")
    ).dt.total_seconds() / 86400.0

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["longitude"], df["latitude"]),
        crs="EPSG:4326",
    )

    # Save outputs
    gpkg_path = output_dir / "ridgecrest_catalog.gpkg"
    csv_path = output_dir / "ridgecrest_catalog.csv"
    gdf.to_file(gpkg_path, driver="GPKG")
    df.to_csv(csv_path, index=False)
    log.info(f"Saved catalog → {gpkg_path}")

    _log_catalog_summary(gdf)
    return gdf


def download_fault_traces(output_dir: Path = Path("data")) -> gpd.GeoDataFrame:
    """
    Download USGS Quaternary Fault and Fold Database fault traces.

    Returns
    -------
    GeoDataFrame
        Fault traces clipped to the Ridgecrest study region.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Downloading USGS Quaternary Fault Database...")
    faults_path = output_dir / "qfaults.zip"

    response = requests.get(QFFDB_URL, timeout=120, stream=True)
    response.raise_for_status()
    with open(faults_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    faults = gpd.read_file(f"zip://{faults_path}")
    faults = faults.to_crs("EPSG:4326")

    # Clip to study region
    bbox = (
        STUDY_REGION["minlongitude"],
        STUDY_REGION["minlatitude"],
        STUDY_REGION["maxlongitude"],
        STUDY_REGION["maxlatitude"],
    )
    faults_clip = faults.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]

    out_path = output_dir / "ridgecrest_faults.gpkg"
    faults_clip.to_file(out_path, driver="GPKG")
    log.info(f"Saved {len(faults_clip)} fault features → {out_path}")
    return faults_clip


def _log_catalog_summary(gdf: gpd.GeoDataFrame) -> None:
    log.info("─── Catalog Summary ───────────────────────────────")
    log.info(f"  Total events:        {len(gdf):,}")
    log.info(f"  Magnitude range:     {gdf['magnitude'].min():.1f} – {gdf['magnitude'].max():.1f}")
    log.info(f"  Depth range (km):    {gdf['depth_km'].min():.1f} – {gdf['depth_km'].max():.1f}")
    log.info(f"  Time span:           {gdf['time'].min()} → {gdf['time'].max()}")
    log.info(f"  M≥5.0 events:        {(gdf['magnitude'] >= 5.0).sum()}")
    log.info(f"  M≥3.0 events:        {(gdf['magnitude'] >= 3.0).sum()}")
    log.info("────────────────────────────────────────────────────")


if __name__ == "__main__":
    gdf = download_catalog()
    download_fault_traces()
