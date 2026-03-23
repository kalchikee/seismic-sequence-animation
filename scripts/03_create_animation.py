#!/usr/bin/env python3
"""
03_create_animation.py
----------------------
Generate the MP4 animation of seismicity evolution and export key frames.

Reads: data/ridgecrest_catalog_analyzed.gpkg, data/ridgecrest_faults.gpkg
Writes: results/ridgecrest_sequence.mp4, results/frames/*.png

Requirements: ffmpeg must be installed
  conda: conda install ffmpeg
  brew:  brew install ffmpeg
  apt:   sudo apt install ffmpeg
"""

import logging
import shutil
import sys
from pathlib import Path

import geopandas as gpd

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.animation import create_animation, export_key_frames

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def main() -> None:
    # Check for ffmpeg
    if shutil.which("ffmpeg") is None:
        log.error(
            "ffmpeg not found. Install it first:\n"
            "  conda: conda install ffmpeg\n"
            "  brew:  brew install ffmpeg\n"
            "  apt:   sudo apt install ffmpeg"
        )
        sys.exit(1)

    data_dir = Path("data")
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    # Load analyzed catalog
    catalog_path = data_dir / "ridgecrest_catalog_analyzed.gpkg"
    if not catalog_path.exists():
        log.error(f"Run 02_analyze_sequence.py first: {catalog_path}")
        sys.exit(1)

    log.info("Loading catalog...")
    catalog = gpd.read_file(catalog_path)

    # Load faults if available
    fault_path = data_dir / "ridgecrest_faults.gpkg"
    faults = gpd.read_file(fault_path) if fault_path.exists() else None

    # Create main animation
    log.info("Creating animation (this may take a few minutes)...")
    create_animation(
        catalog=catalog,
        faults=faults,
        output_path=results_dir / "ridgecrest_sequence.mp4",
        frames_per_day=1,
    )

    # Export key frames
    log.info("Exporting key frames...")
    export_key_frames(
        catalog=catalog,
        faults=faults,
        output_dir=results_dir / "frames",
        days=[1, 7, 30, 90, 184],
    )

    print("\n✓ Animation complete:")
    print(f"  Video → {results_dir}/ridgecrest_sequence.mp4")
    print(f"  Frames → {results_dir}/frames/")


if __name__ == "__main__":
    main()
