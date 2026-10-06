# Data dictionary

Grid: ISEA3H resolution 16 (DGGRID), WGS84 (EPSG:4326). Areas in km², population in persons.

**Spatial domain (`S60_N90`):** every file covers the land cells whose centroid lies at or north
of 60°S — 114,791,301 cells. Antarctica and other land south of 60°S are not included.

**Gzip CSV conventions** (`*.csv.gz`): one row per cell, sorted by `seqnum`; empty fields are
missing values; pixel counts and `gaul0_code` are integers; flags are 1 (true) or 0 (false).
Decimals written: 7 for coordinates and areas (about 1 cm and 0.1 m²), 5 for fractions (one 10 m
pixel is ~0.000085 of a cell's land), 3 for distances, population, densities and `n_valid_pop_px`.
Very small values may be written in scientific notation (e.g. `8.5e-05`).

## Grid

### `dgg_isea3h16_land_shapefiles_S60_N90.zip`

377 zipped shapefiles, one per layer: `dgg_isea3h16_land_<n|s>LL_<e|w>LLL.zip` (10° × 10° tile
named by its south-west corner). n = 114,791,301.

| Field | Type | Description |
|---|---|---|
| `seqnum` | text | cell id (DGGRID sequence number), stored as text because 9-digit integers overflow shapefile numeric fields |
| geometry | polygon / multipolygon | cell boundary; cells crossing ±180° are split into east and west parts |

## `dgg_isea3h16_land_centroids_iso3_S60_N90.csv.gz` — cell centroids and country assignment

Gzip CSV, one row per grid cell (114,791,301), sorted by `seqnum`.

| Column | Description |
|---|---|
| `seqnum` | cell id |
| `longitude`, `latitude` | cell centroid (degrees) |
| `gaul0_code`, `gaul0_name`, `iso3_code` | GAUL 2024 level 0 unit matched to the cell (empty if none) |
| `iso3_admin` | de facto country used in the country summaries (empty = unassigned) |
| `admin_source` | how `iso3_admin` was set: `gaul` (GAUL code as published), `crosswalk` (disputed GAUL feature, de facto administrator), `naturalearth` (Jammu & Kashmir / India–China border areas, per cell from Natural Earth), `naturalearth_remote` (> 50 km from GAUL, from Natural Earth), `unassigned` (Abyei, Bir Tawil, Spratly Islands), `beyond_max_dist` (> 50 km from GAUL and not in Natural Earth) |
| `match_type` | `intersects` (centroid inside a GAUL polygon), `nearest` (nearest GAUL polygon), `ne_remote` (Natural Earth, remote islands) |
| `dist_km` | geodesic distance from the centroid to the matched polygon for `nearest` / `ne_remote` matches (km) |
| `wrapped` | 1 if the cell crosses the antimeridian (stored as east/west parts) |

## `dgg_cells_shared_landscapes_S60_N90_dataset.csv.gz` — per-cell results

Gzip CSV, one row per cell with centroid latitude ≥ 60°S (114,791,301), sorted by `seqnum`;
29 columns. `popland` and `shland` were computed from unrounded values, so a cell whose fraction
lies within 0.000005 of a threshold may appear to contradict its flag.

**Location and country**

| Column | Description |
|---|---|
| `seqnum` | cell id |
| `longitude`, `latitude` | cell centroid (degrees) |
| `gaul0_code`, `gaul0_name`, `iso3_code`, `iso3_admin` | as in `dgg_isea3h16_land_centroids_iso3_S60_N90.csv.gz` |

**Earth Engine pixel counts** — numbers of 10 m pixels (EPSG:4326, 10 m scale) whose centres fall
inside the cell; definitions in `docs/METHODS.md`.

| Column | Description |
|---|---|
| `n_total` | valid pixels: not permanent water in ESA WorldCover 2021 and ≥ 1 Dynamic World observation in 2024 |
| `n_land` | valid pixels that are not persistent freshwater or persistent snow/ice |
| `n_freshwater` | persistent freshwater pixels (Dynamic World water ≥ 95% of observations, clusters ≥ 10 pixels) |
| `n_ice` | persistent snow/ice pixels (Dynamic World snow/ice ≥ 99% of observations, clusters ≥ 10 pixels) |
| `n_crop`, `n_built`, `n_barren` | land pixels where that Dynamic World class occurs in > 40% of observations |
| `n_pasture` | land pixels of cultivated grassland (Global Pasture Watch, 2022) |
| `n_other_intensive_treecrop` | land pixels of planted trees, class 2 (Spatial Database of Planted Trees, Version 2.0) |
| `n_nonhabitat` | land pixels that are crop, built-up, pasture or tree crop |
| `n_habitat` | `n_land − n_nonhabitat` |

**Population (WorldPop R2025A, 2024, 100 m constrained)**

| Column | Description |
|---|---|
| `cell_area_km2` | cell area (equal-area projection EPSG:6933) |
| `pop_sum` | persons in the cell: sum of pixel values × fraction of each pixel covered by the cell (empty outside the raster's latitude extent) |
| `pop_density` | `pop_sum / cell_area_km2` |
| `n_valid_pop_px` | number of populated 100 m pixels covering the cell (coverage-weighted, so not necessarily a whole number) |

**Areas and classification**

| Column | Description |
|---|---|
| `land_area_km2` | land area: `n_land` × area of one 10 m pixel at the cell's latitude (WGS84 ellipsoid) |
| `habitat_area_km2` | habitat area, likewise from `n_habitat` |
| `nonhabitat_frac` | non-habitat fraction of the cell's land, `n_nonhabitat / n_land` (fraction, 0-1) |
| `habitat_frac` | habitat fraction of the cell's land, `n_habitat / n_land` (= `habitat_area_km2 / land_area_km2`; fraction, 0-1) |
| `pop_density_land` | persons per km² of land, `pop_sum / land_area_km2` (used for the classification) |
| `popland` | Populated Landscape (1/0): `n_land` > 0, `nonhabitat_frac` > 0.01 and `pop_density_land` ≥ 1 |
| `shland` | Shared Landscape (1/0): `popland` and `habitat_frac` ≥ 0.20 |

## `dataset_summaries.zip` — summary tables

Contains `global_summary.csv` (one row, `group` = `Global`) and `country_summary.csv` (one row per
de facto country: `group` = name, `iso3`; `--` = unassigned), which share these columns. Population fractions use
all cells with centroid latitude ≥ 60°S (WorldPop R2025A, 2024); fractions are rounded to 4 decimals.

| Column | Description |
|---|---|
| `land_area_km2` | land (area, km2)|
| `popland_area_km2` | land in Populated Landscapes (area, km2)|
| `total_habitat_area_km2` | habitat (area, km2)|
| `habitat_in_popland_area_km2` | habitat in Populated Landscapes (area, km2)|
| `shland_area_km2` | land in Shared Landscapes (area, km2)|
| `shland_habitat_area_km2` | habitat in Shared Landscapes (area, km2)|
| `n_hex` | cells with land (`n_land` > 0) |
| `frac_land_popland` | fraction of land that is Populated Landscapes (fraction, 0-1)|
| `frac_habitat_in_popland` | fraction of habitat in Populated Landscapes (fraction, 0-1)|
| `frac_shland` | fraction of land that is Shared Landscapes (fraction, 0-1)|
| `frac_shland_in_popland` | fraction of Populated Landscapes that is Shared Landscapes (fraction, 0-1)|
| `frac_population_in_shland` | fraction of total population residing in Shared Landscapes (fraction, 0-1)|
