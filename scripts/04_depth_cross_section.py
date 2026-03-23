#!/usr/bin/env python3
"""
04_depth_cross_section.py
-------------------------
Build fault-perpendicular depth cross-section of the Ridgecrest aftershock
sequence, revealing seismogenic zone structure and brittle-ductile transition.
"""

import logging
import sys
from pathlib import Path

import geopandas as gpd

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.cross_section import plot_cross_section

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


def main() -> None:
    gpkg = DATA_DIR / "ridgecrest_catalog_analyzed.gpkg"
    if not gpkg.exists():
        print("Run 02_analyze_sequence.py first.")
        sys.exit(1)

    gdf = gpd.read_file(gpkg)
    plot_cross_section(gdf, output_path=RESULTS_DIR / "depth_cross_section.png")
    print(f"✓ Cross-section → {RESULTS_DIR}/depth_cross_section.png")


if __name__ == "__main__":
    main()
