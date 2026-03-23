#!/usr/bin/env python3
"""
01_download_catalog.py
----------------------
Entry-point script: download the 2019 Ridgecrest earthquake catalog
from the USGS FDSN web service and save to data/.

Usage:
    python scripts/01_download_catalog.py
    python scripts/01_download_catalog.py --minmagnitude 1.0 --endtime 2020-06-01
"""

import argparse
import logging
import sys
from pathlib import Path

# Make src importable from scripts/
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.catalog_download import download_catalog, download_fault_traces

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download Ridgecrest earthquake catalog from USGS ComCat"
    )
    parser.add_argument("--starttime", default="2019-07-01",
                        help="Start date (ISO 8601, default: 2019-07-01)")
    parser.add_argument("--endtime", default="2020-01-01",
                        help="End date (ISO 8601, default: 2020-01-01)")
    parser.add_argument("--minmagnitude", type=float, default=1.5,
                        help="Minimum magnitude (default: 1.5)")
    parser.add_argument("--output-dir", default="data",
                        help="Output directory (default: data/)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)

    print("=" * 60)
    print("  Ridgecrest Sequence Catalog Download")
    print("=" * 60)
    print(f"  Time window : {args.starttime} → {args.endtime}")
    print(f"  Min magnitude: M{args.minmagnitude}")
    print(f"  Output dir  : {output_dir.resolve()}")
    print()

    # Download earthquake catalog
    catalog = download_catalog(
        starttime=args.starttime,
        endtime=args.endtime,
        minmagnitude=args.minmagnitude,
        output_dir=output_dir,
    )

    # Download fault traces
    print("\nDownloading USGS Quaternary Fault traces...")
    faults = download_fault_traces(output_dir=output_dir)

    print(f"\n✓ Downloaded {len(catalog):,} earthquakes")
    print(f"✓ Downloaded {len(faults):,} fault features")
    print(f"✓ Data saved to: {output_dir.resolve()}/")
    print("\nNext step: python scripts/02_analyze_sequence.py")


if __name__ == "__main__":
    main()
