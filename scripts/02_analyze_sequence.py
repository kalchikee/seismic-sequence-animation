#!/usr/bin/env python3
"""
02_analyze_sequence.py
----------------------
Run spatiotemporal analysis on the Ridgecrest catalog:
  - Compute along-strike/along-dip distances
  - Fit Omori decay law
  - Compute completeness magnitude and b-values
  - Analyze spatial migration patterns

Reads: data/ridgecrest_catalog.gpkg
Writes: results/ (CSV files with analysis outputs)
"""

import logging
import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.spatiotemporal_analysis import (
    compute_fault_distances,
    fit_omori_law,
    completeness_magnitude,
    aki_bvalue,
    spatial_bvalue_map,
    migration_analysis,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


def main() -> None:
    # ── Load catalog ───────────────────────────────────────────────────────────
    gpkg = DATA_DIR / "ridgecrest_catalog.gpkg"
    if not gpkg.exists():
        log.error(f"Catalog not found: {gpkg}. Run 01_download_catalog.py first.")
        sys.exit(1)

    log.info(f"Loading catalog from {gpkg}")
    gdf = gpd.read_file(gpkg)
    log.info(f"  {len(gdf):,} events loaded")

    # ── Fault distances ────────────────────────────────────────────────────────
    log.info("Computing fault plane distances...")
    gdf = compute_fault_distances(gdf)

    # ── Completeness magnitude ─────────────────────────────────────────────────
    mc = completeness_magnitude(gdf["magnitude"].values)
    log.info(f"Completeness magnitude Mc = {mc}")

    # ── Overall b-value ────────────────────────────────────────────────────────
    b, sigma_b = aki_bvalue(gdf["magnitude"].values, mc)
    log.info(f"Overall b-value: {b:.3f} ± {sigma_b:.3f} (Mc={mc})")

    # ── Omori law ──────────────────────────────────────────────────────────────
    log.info("Fitting Omori decay law...")
    omori = fit_omori_law(gdf)
    omori_df = pd.DataFrame({
        "t_days": omori["t_bins"],
        "rate_observed": omori["rate_observed"],
        "rate_fitted": omori["rate_fitted"],
    })
    omori_df.to_csv(RESULTS_DIR / "omori_fit.csv", index=False)

    params_df = pd.DataFrame([{
        "K": omori["K"], "K_std": omori["K_std"],
        "c": omori["c"], "c_std": omori["c_std"],
        "p": omori["p"], "p_std": omori["p_std"],
        "Mc": mc, "b": b, "sigma_b": sigma_b,
    }])
    params_df.to_csv(RESULTS_DIR / "sequence_parameters.csv", index=False)
    log.info(f"  Omori parameters → {RESULTS_DIR}/sequence_parameters.csv")

    # ── Spatial b-value map ────────────────────────────────────────────────────
    log.info("Computing spatial b-value map...")
    bval_map = spatial_bvalue_map(gdf, bin_deg=0.1, min_events=30)
    bval_map.to_csv(RESULTS_DIR / "bvalue_spatial.csv", index=False)
    log.info(f"  {len(bval_map)} spatial bins → {RESULTS_DIR}/bvalue_spatial.csv")

    # ── Migration analysis ─────────────────────────────────────────────────────
    log.info("Analyzing seismicity migration...")
    migration = migration_analysis(gdf)
    mig_df = pd.DataFrame([migration])
    mig_df.to_csv(RESULTS_DIR / "migration_analysis.csv", index=False)

    # ── Save enriched catalog ──────────────────────────────────────────────────
    gdf.to_file(DATA_DIR / "ridgecrest_catalog_analyzed.gpkg", driver="GPKG")
    log.info("Enriched catalog saved → data/ridgecrest_catalog_analyzed.gpkg")

    # ── Summary ────────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Analysis Results Summary")
    print("=" * 60)
    print(f"  Total events (M≥Mc):  {(gdf['magnitude'] >= mc).sum():,}")
    print(f"  Completeness Mc:      {mc}")
    print(f"  b-value:              {b:.3f} ± {sigma_b:.3f}")
    print(f"  Omori p:              {omori['p']:.3f} ± {omori['p_std']:.3f}")
    print(f"  Omori K:              {omori['K']:.1f} ± {omori['K_std']:.1f}")
    print(f"  Migration model:      {migration['preferred_model']}")
    print(f"  Migration velocity:   {migration['velocity_km_per_day']:.2f} km/day")
    print()
    print("  Next step: python scripts/03_create_animation.py")


if __name__ == "__main__":
    main()
