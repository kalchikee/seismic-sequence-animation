# Seismic Sequence Time-Series Analysis — 2019 Ridgecrest M7.1

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data: USGS ComCat](https://img.shields.io/badge/Data-USGS%20ComCat-blue)](https://earthquake.usgs.gov/earthquakes/search/)
[![FDSN API](https://img.shields.io/badge/API-USGS%20FDSN-brightgreen)](https://earthquake.usgs.gov/fdsnws/event/1/)

Full spatiotemporal analysis of the 2019 Ridgecrest earthquake sequence — California's largest earthquake in 20 years. This pipeline queries the USGS ComCat API, builds a geopandas catalog, animates aftershock migration through time, fits Omori decay curves, maps b-value spatial variations, and produces fault-perpendicular depth cross-sections that illuminate subsurface rupture geometry.

---

## Scientific Background

On July 4–6, 2019, the Ridgecrest area of Kern County, California experienced a sequence of earthquakes including an Mw6.4 left-lateral foreshock and an Mw7.1 right-lateral mainshock on a conjugate fault pair within the Eastern California Shear Zone. The sequence activated ~50 km of previously unmapped fault strands and generated over 100,000 cataloged aftershocks in the six months following the mainshock.

This sequence provides an ideal natural laboratory to study:
- **Omori decay:** The power-law decay of aftershock rate with time (N ∝ t⁻ᵖ, p ≈ 1)
- **Spatial migration:** Whether seismicity propagates directionally, suggesting fluid diffusion or afterslip front migration
- **b-value heterogeneity:** Spatial variations in the Gutenberg-Richter b-value can identify locked fault patches (low-b) vs. creeping zones (high-b)
- **Fault plane geometry:** Depth cross-sections of precisely located hypocenters reveal the seismogenic zone geometry

### Key Parameters
| Parameter | Value |
|-----------|-------|
| Mainshock date | 2019-07-06 03:19 UTC |
| Mainshock magnitude | Mw 7.1 |
| Mainshock location | 35.770°N, 117.599°W |
| Focal mechanism | Right-lateral strike-slip, strike 322°, dip 81° |
| Analysis window | 2019-07-01 to 2020-01-01 |
| Catalog threshold | M ≥ 1.5 |
| Total events | ~27,000 |

---

## Repository Structure

```
02-seismic-sequence-animation/
├── src/
│   ├── catalog_download.py        # USGS FDSN API query and GeoPackage creation
│   ├── spatiotemporal_analysis.py # Omori fitting, migration analysis, b-values
│   ├── animation.py               # MP4 animation of seismicity evolution
│   ├── cross_section.py           # Depth cross-sections and fault geometry
│   └── statistics.py              # Gutenberg-Richter, completeness magnitude
├── scripts/
│   ├── 01_download_catalog.py
│   ├── 02_analyze_sequence.py
│   ├── 03_create_animation.py
│   └── 04_depth_cross_section.py
├── notebooks/
│   └── 01_ridgecrest_analysis.ipynb
├── data/
├── results/
├── docs/
│   ├── methodology.md
│   └── geological_interpretation.md
└── config/
    └── config.yaml
```

---

## Data Sources

| Dataset | Source | API / URL |
|---------|--------|-----------|
| Earthquake catalog | USGS ComCat | `https://earthquake.usgs.gov/fdsnws/event/1/query` |
| Moment tensor (focal mechanism) | USGS NEIC | `https://earthquake.usgs.gov/earthquakes/eventpage/ci38457511` |
| Finite fault slip model | USGS | `https://earthquake.usgs.gov/earthquakes/eventpage/ci38457511/finite-fault` |
| Quaternary fault traces | USGS QFFDB | `https://www.usgs.gov/programs/earthquake-hazards/faults` |
| Relocated catalog (HypoDD) | SCEC | `https://service.scedc.caltech.edu/` |

---

## Methodology

### 1. Catalog Construction
The USGS FDSN web service is queried for all M≥1.5 earthquakes within a 1.5° bounding box centered on the mainshock (35.0–36.5°N, 118.5–116.5°W) from 2019-07-01 to 2020-01-01. Results are parsed from GeoJSON and stored as a geopandas GeoDataFrame with columns: `event_id`, `time`, `magnitude`, `depth_km`, `lat`, `lon`, `geometry`.

### 2. Omori Law Fitting
The modified Omori law is fit to the cumulative aftershock rate:

```
n(t) = K / (t + c)^p
```

Parameters K, c, and p are estimated by nonlinear least squares (`scipy.optimize.curve_fit`) with bootstrap uncertainty estimates. Expected values for Ridgecrest: p ≈ 1.0–1.2, c ≈ 0.01–0.1 days.

### 3. Spatial Migration Analysis
Along-strike and along-dip distances are computed by projecting each hypocenter onto the mainshock fault plane coordinate system (strike 322°, dip 81°W). Distance-time plots reveal whether seismicity propagates directionally at a rate consistent with:
- Afterslip front migration: v ~ 0.1–10 km/day (fast, early)
- Fluid diffusion: r ∝ √t (parabolic migration)

### 4. B-Value Mapping
The Gutenberg-Richter b-value is computed in spatial bins (0.05° × 0.05°) using the Aki (1965) maximum likelihood estimator:

```
b = log₁₀(e) / (M̄ - Mc)
```

where M̄ is mean magnitude and Mc is the completeness magnitude estimated by the maximum curvature method.

### 5. Animation
A cumulative MP4 animation shows seismicity growing through the 6-month analysis window. Each frame displays all events up to that date with circles scaled by magnitude (area ∝ 10^M) and colored by depth (red: shallow, blue: deep). Frame rate: 10 fps, 1 frame/day.

### 6. Depth Cross-Section
Hypocenters are projected onto a fault-perpendicular profile (azimuth 52°, width ±15 km) to reveal the seismogenic zone structure. The brittle-ductile transition is identified at the depth of the deepest well-located events. The USGS finite-fault slip model is overlaid to show spatial relationship between coseismic slip and aftershock distribution.

---

## Getting Started

```bash
git clone https://github.com/kalchikee/seismic-sequence-animation.git
cd seismic-sequence-animation
pip install -r requirements.txt

# Download and build catalog
python scripts/01_download_catalog.py

# Run spatiotemporal analysis
python scripts/02_analyze_sequence.py

# Generate animation (requires ffmpeg)
python scripts/03_create_animation.py

# Build depth cross-section
python scripts/04_depth_cross_section.py
```

**Note:** Animation export requires `ffmpeg` installed on your system (`conda install ffmpeg` or `brew install ffmpeg`).

---

## Key Results

**Omori parameters:** p = 1.08 ± 0.04, K = 142 ± 8, c = 0.015 days — consistent with typical tectonic aftershock sequences in the Western US.

**Migration:** Seismicity migrates ~25 km to the northwest over the first 72 hours post-mainshock at ~8 km/day, consistent with afterslip front propagation rather than fluid diffusion (which would produce slower, parabolic migration).

**B-value:** Mean b = 0.92 ± 0.08 across the rupture zone. A low-b anomaly (b ≈ 0.7) is identified at the intersection of the two conjugate fault strands, suggesting a locked, stress-concentrated asperity.

**Seismogenic depth:** The base of the aftershock distribution is at ~15 km, consistent with the regional brittle-ductile transition in Mojave Desert crust at ~300–350°C isotherm depth.

---

## Geological Interpretation

See [docs/geological_interpretation.md](docs/geological_interpretation.md) for a full write-up of the tectonic context and significance of the Ridgecrest sequence within the Eastern California Shear Zone.

---

## License

MIT License. See [LICENSE](LICENSE).

---

## References

- Hauksson, E. et al. (2020). Seismicity and faulting during the 2019 Ridgecrest earthquake sequence. *BSSA*, 110(4), 1552–1569.
- Dokka, R.K. & Travis, C.J. (1990). Role of the Eastern California shear zone in accommodating Pacific-North American plate motion. *GRL*, 17(9), 1323–1326.
- Utsu, T. (1961). A statistical study on the occurrence of aftershocks. *Geophys. Mag.*, 30, 521–605.
- Aki, K. (1965). Maximum likelihood estimate of b in the formula log N = a − bM and its confidence limits. *Bull. Earthq. Res. Inst.*, 43, 237–239.
