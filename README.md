# Habitat Areas in Populated Landscapes — A Global DGG Computational Pipeline

Code to reproduce the data for 'A global dataset of potential habitat areas in populated landscapes' on a global ~1 km² hexagonal grid:

1. build a **global all-land discrete global grid** (ISEA3H resolution 16, 1.18491 km² cells);
2. count **10 m land-cover pixels** per cell in Google Earth Engine (Dynamic World 2024, ESA
   WorldCover 2021, Global Pasture Watch, Spatial Database of Planted Trees);
3. add **country** (GAUL 2024, de facto policy for disputed areas) and **population**
   (WorldPop R2025A 2024, 100 m) to every cell;
4. calculate habitat and landscape extents (Populated / Shared) and write **global and country summaries**.

Data: **Zenodo DOI: [10.5281/zenodo.23168754](https://doi.org/10.5281/zenodo.23168754)** — grid, GEE counts, per-cell
results and summary tables (see [Data on Zenodo](#data-on-zenodo)).

---

## Pipeline

| Step | Script | Runs in | Output |
|---|---|---|---|
| 0 | `R/00_install_packages.R` | R | R packages |
| 1 | `R/01_build_dgg.R` | R (~5 h) | `data/dgg/` — global land DGG, one layer per 10° tile |
| 2 | `R/02_validate_dgg.R` | R (minutes) | validation report (checks 0–7) |
| 3 | `gee/03_extract_counts.js`, `python/03b_combine_counts.py` | Earth Engine batch (377 export tasks; ~1.1 million EECU-hours); Python | pixel-count CSVs per layer, combined into `data/gee_counts/dgg_counts_2024.csv.gz` |
| 4 | `python/00_prepare_gaul_l0.py`, `python/04_download_worldpop.py` | Python | GAUL 2024 level 0 GeoPackage; WorldPop 100 m country rasters + global VRT |
| 5 | `python/05_postprocess.py` | Python (~4–5 h first run; cached afterwards) | per-cell tables + summaries in `data/outputs/` |

All Python settings are in **`python/config.py`**; R settings are at the top of each R script.
Run everything **from the repository root**. Set `DATA_DIR` to keep `data/` elsewhere
(e.g. a fast local disk; avoid synced folders such as Dropbox for the multi-GB intermediates).

**Rerunning steps 4–5 without Earth Engine:** build the grid (steps 1–2), then set `GEE_INPUT`
in `python/config.py` to the published `dgg_cells_shared_landscapes_S60_N90_dataset.csv.gz`, which
contains the Earth Engine counts.

---

## Setup

```bash
conda env create -f environment.yml
conda activate populated-landscapes-dgg
Rscript R/00_install_packages.R
```

Tested on macOS (Apple silicon, 128 GB RAM): GEOS 3.13.0, GDAL 3.8.5, PROJ 9.5.1, Python 3.10.
For exact reproducibility pin your environment (`conda env export > environment.lock.yml`;
`renv::snapshot()` for R) and commit the lock files.

## Inputs (place in `data/inputs/`)

| Input | Used by | Source |
|---|---|---|
| Natural Earth 10m land (`ne_10m_land/`) | steps 1–2 | [Natural Earth physical vectors](https://www.naturalearthdata.com/downloads/10m-physical-vectors/) (public domain) |
| OSM land polygons, WGS84 split (`land-polygons-split-4326/`) | steps 1–2 | [osmdata.openstreetmap.de](https://osmdata.openstreetmap.de/data/land-polygons.html) (ODbL, © OpenStreetMap contributors) — record the "Last update" date you use; the reference build used 2026-09-30T03:38 |
| GAUL 2024 level 1 (`GAUL_2024_L1/`), dissolved to level 0 (`gaul2024_l0.gpkg`) by `python/00_prepare_gaul_l0.py` | step 5 | FAO Global Administrative Unit Layers 2024, [level 1](https://data.apps.fao.org/catalog/iso/34f97afc-6218-459a-971d-5af1162d318a) (CC BY 4.0 and the [GAUL 2024 Terms of Use](https://data.apps.fao.org/catalog/dataset/77d2cb0d-348a-4fdf-aeec-0f47c5d14f28/resource/a9067338-7bcd-4d82-8b65-adb1d0277eb3/download/gaul2024termsofuse.pdf)); level 0 by dissolving on `gaul0_code` as in the [GEE community catalog tutorial](https://gee-community-catalog.org/tutorials/examples/gaul_aggregate/#section-3-handling-gee-limitations) (section 4). Cite as: FAO. 2024. Global Administrative Unit Layers (GAUL). [Accessed on 8 July 2026]. https://data.apps.fao.org/?lang=en. Licence: CC-BY-4.0 |
| Natural Earth 10m Admin 0 countries (`ne_10m_admin_0_countries/`) | step 5 | [Natural Earth cultural vectors](https://www.naturalearthdata.com/downloads/10m-cultural-vectors/) (public domain) |
| WorldPop R2025A 100 m constrained, 2024 | step 5 | downloaded by step 4 from [WorldPop](https://hub.worldpop.org/geodata/listing?id=135) |

Earth Engine datasets (step 3):

| Dataset | Earth Engine ID | Used for | Citation | Licence / terms |
|---|---|---|---|---|
| Dynamic World V1 (2024) | `GOOGLE/DYNAMICWORLD/V1` | valid observations; water, snow/ice, crop, built, bare | Brown, C.F., Brumby, S.P., Guzder-Williams, B. et al. Dynamic World, Near real-time global 10 m land use land cover mapping. *Sci Data* 9, 251 (2022). [doi:10.1038/s41597-022-01307-4](https://doi.org/10.1038/s41597-022-01307-4) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Required attribution: "This dataset is produced for the Dynamic World Project by Google in partnership with National Geographic Society and the World Resources Institute." Contains modified Copernicus Sentinel data [2015-present]; see the [Sentinel Data Legal Notice](https://sentinels.copernicus.eu/documents/247904/690755/Sentinel_Data_Legal_Notice). |
| ESA WorldCover 10 m 2021 v200 | `ESA/WorldCover/v200/2021` | permanent-water (ocean) mask | Zanaga, D., Van De Kerchove, R., Daems, D., De Keersmaecker, W., Brockmann, C., Kirches, G., Wevers, J., Cartus, O., Santoro, M., Fritz, S., Lesiv, M., Herold, M., Tsendbazar, N.E., Xu, P., Ramoino, F., Arino, O. (2022). ESA WorldCover 10 m 2021 v200. [doi:10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221) | [CC BY 4.0](https://spdx.org/licenses/CC-BY-4.0.html) |
| Global Pasture Watch, grassland class v1 (2022) | `projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c` | cultivated grassland (pasture) | Parente, L., Sloat, L., Mesquita, V., et al. (2024). Global Pasture Watch - Annual grassland class and extent maps at 30-m spatial resolution (2000—2022) (Version v1) [Data set]. Zenodo. [doi:10.5281/zenodo.13890401](https://doi.org/10.5281/zenodo.13890401); Parente, L., Sloat, L., Mesquita, V., et al. (2024). Annual 30-m maps of global grassland class and extent (2000–2022) based on spatiotemporal Machine Learning. *Scientific Data*. [doi:10.1038/s41597-024-04139-6](http://doi.org/10.1038/s41597-024-04139-6) | [CC BY 4.0](https://spdx.org/licenses/CC-BY-4.0.html) |
| Spatial Database of Planted Trees (SDPT) Version 2.0 | `projects/ee-duzhenrong02/assets/SDPTV3` (community-hosted asset containing SDPT Version 2.0) | tree crops (class 2) | Richter, J., E. Goldman, N. Harris, D. Gibbs, M. Rose, S. Peyer, S. Richardson, and H. Velappan. 2024. "Spatial Database of Planted Trees (SDPT Version 2.0)." Technical Note. Washington, DC: World Resources Institute. [doi:10.46830/writn.23.00073](https://doi.org/10.46830/writn.23.00073) | [CC BY 4.0](https://spdx.org/licenses/CC-BY-4.0.html) |

---

## Running

### 1–2. Grid

```bash
Rscript R/01_build_dgg.R  2>&1 | tee build_dgg.log
Rscript R/02_validate_dgg.R 2>&1 | tee validate_dgg.log
```

Every check must report `OK`. The grid is **global** (Antarctica included); the study
domain (centroid latitude ≥ 60°S) is applied in steps 3 and 5.

### 3. Earth Engine counts

1. Upload each layer zip from `data/dgg/` whose tile latitude is ≥ −60 (`n00…n80`, `s10…s60`,
   and `north_cap` if present; not `s70…s90` or `south_cap`) as a **table asset** into one
   asset folder.
   Many files: stage them in Cloud Storage and use `earthengine upload table` (see the end of
   `R/01_build_dgg.R`).
2. Set `ASSET_FOLDER` ([PLACEHOLDER] in `gee/03_extract_counts.js`), run it, and start the export tasks. Process
   the layer list in slices (`FIRST_LAYER`, `N_LAYERS`) and confirm every task completed.
3. Download all CSVs from the Drive folder into `data/gee_counts/raw/`, then combine them:

```bash
python -u python/03b_combine_counts.py
```

This writes `data/gee_counts/dgg_counts_2024.csv.gz` (one row per cell, integer counts), the table
read by step 5. Alternatively, set `GEE_INPUT` in `python/config.py` to a folder to read every
file matching `GEE_INPUT_PATTERN` directly.

### 4–5. Population, countries, classification

```bash
python -u python/00_prepare_gaul_l0.py
python -u python/04_download_worldpop.py 2>&1 | tee download_worldpop.log
python -u python/05_postprocess.py      2>&1 | tee postprocess.log
```

`05_postprocess.py` caches its expensive steps in `data/work/` and resumes after interruption.
It must run from a terminal (it starts worker processes), not a notebook.

---

## Outputs (`data/outputs/`)

- `dgg_isea3h16_land_centroids_iso3_S60_N90.csv.gz` — centroid and country assignment of every
  cell with centroid latitude ≥ 60°S (GAUL fields, de facto `iso3_admin`, how it was assigned).
- `dgg_cells_shared_landscapes_S60_N90_dataset.csv.gz` — per-cell results (29 columns): location,
  country, Earth Engine pixel counts, population, land and habitat areas and fractions,
  `popland` (Populated Landscape), `shland` (Shared Landscape).
- `dataset_summaries.zip` — `global_summary.csv` (one row) and `country_summary.csv`: areas,
  fractions and the fraction of population in Shared Landscapes (also in `dataset_summaries/`).
- `dgg_cells_country_population.parquet` — intermediate table of step 5.

Reference run (OSM 2026-09-30, Dynamic World 2024, WorldPop 2024): 114,791,301 domain cells;
land 127.61 million km²; Populated Landscapes 26.14 million km² (20.5% of land); Shared
Landscapes 18.87 million km² (14.8% of land); 100.00% of WorldPop population captured.

## Method summary

See [`docs/METHODS.md`](docs/METHODS.md).

## Data on Zenodo

| File | Content |
|---|---|
| `dgg_isea3h16_land_shapefiles_S90_N90.zip` | global land DGG: 474 zipped layer shapefiles (`seqnum` as text), 126,509,398 cells |
| `dgg_isea3h16_land_centroids_iso3_S60_N90.csv.gz` | centroid and country assignment of each cell with centroid latitude ≥ 60°S |
| `dgg_cells_shared_landscapes_S60_N90_dataset.csv.gz` | per-cell results: Earth Engine counts, population, land and habitat areas, landscape classes |
| `dataset_summaries.zip` | `global_summary.csv` and `country_summary.csv` |
| `DATA_DICTIONARY.md` | column definitions for all files (also in [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md)) |

`S90_N90` marks the global grid (all latitudes); `S60_N90` marks files covering cells with
centroid latitude ≥ 60°S. The shapefile archive bundles the layer zips written by step 1:
`cd data/dgg && zip -0 dgg_isea3h16_land_shapefiles_S90_N90.zip dgg_isea3h16_land_*.zip`. Cell
centroids for the full grid are written by step 1 (`*_pts.csv`) or can be computed from the polygons.

WorldPop, GAUL, Natural Earth and OSM inputs are not redistributed; download them from their
sources above.

## Citation and licence

**Data** (Zenodo record): [Open Database License (ODbL) v1.0](https://opendatacommons.org/licenses/odbl/1-0/).
The grid is derived from OpenStreetMap data (© OpenStreetMap contributors, ODbL); the inputs
listed above are credited under their own terms (CC BY 4.0 for GAUL, Dynamic World, ESA
WorldCover, Global Pasture Watch and SDPT; WorldPop and Natural Earth as stated by their
providers).

**Code**: [MIT License](LICENSE).

**Citation**: Ellis, E.C., Depsky, N., Abrams, J.F., Dong, J., Yu, L., Di, Y., Lyu, B., Du, Z.,
Verburg, P., Tapia, H. A global dataset of potential habitat areas in populated landscapes.
*Data in Brief* (2026). [PLACEHOLDER: article DOI, pending]. Data: Zenodo, [doi:10.5281/zenodo.23168754](https://doi.org/10.5281/zenodo.23168754).
Machine-readable metadata: [`CITATION.cff`](CITATION.cff).
