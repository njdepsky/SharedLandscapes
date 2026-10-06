# Data dictionary

Grid: ISEA3H resolution 16 (DGGRID), WGS84 (EPSG:4326). Areas in km², population in persons.
Study domain of the counts and results: cells whose centroid latitude ≥ 60°S (114,791,301 cells).

## Grid

### `dgg_isea3h16_land_v4_shapefiles.zip`

474 zipped shapefiles, one per layer: `dgg_isea3h16_land_<n|s>LL_<e|w>LLL.zip` (10° × 10° tile
named by its south-west corner) and `dgg_isea3h16_land_south_cap.zip`. n = 126,509,398.

| Field | Type | Description |
|---|---|---|
| `seqnum` | text | cell id (DGGRID sequence number), stored as text because 9-digit integers overflow shapefile numeric fields |
| geometry | polygon / multipolygon | cell boundary; cells crossing ±180° are split into east and west parts |

### `dgg_isea3h16_land_v4_centroids.zip`

474 CSV files, `<layer>_pts.csv`, one row per cell.

| Column | Type | Description |
|---|---|---|
| `seqnum` | integer | cell id |
| `longitude`, `latitude` | float | cell centroid, decimal degrees |

### `dgg_isea3h16_land_manifest.csv`

| Column | Description |
|---|---|
| `layer` | layer name |
| `lat0`, `lon0` | south-west corner of the tile (degrees) |
| `n_land_cells` | land cells in the layer (cells of the global all-land grid whose centroid lies in the tile) |
| `seconds` | build time of the layer |

## `dgg_counts_2024.csv.gz` — Earth Engine pixel counts

Gzip CSV, one row per domain cell, sorted by `seqnum16`; all values integers. Counts are numbers of 10 m pixels (EPSG:4326, 10 m scale) whose centres
fall inside the cell; definitions in `docs/METHODS.md`.

| Column | Description |
|---|---|
| `seqnum16` | cell id |
| `n_total` | valid pixels: not permanent water in ESA WorldCover 2021 and ≥ 1 Dynamic World observation in 2024 |
| `n_land` | valid pixels that are not persistent freshwater or persistent snow/ice |
| `n_freshwater` | persistent freshwater pixels (Dynamic World water ≥ 95% of observations, clusters ≥ 10 pixels) |
| `n_ice` | persistent snow/ice pixels (Dynamic World snow/ice ≥ 99% of observations, clusters ≥ 10 pixels) |
| `n_crop`, `n_built`, `n_barren` | land pixels where that Dynamic World class occurs in > 40% of observations |
| `n_pasture` | land pixels of cultivated grassland (Global Pasture Watch, 2022) |
| `n_other_intensive_treecrop` | land pixels of planted trees, class 2 (Spatial Database of Planted Trees, Version 2.0) |
| `n_nonhabitat` | land pixels that are crop, built-up, pasture or tree crop |
| `n_habitat` | `n_land − n_nonhabitat` |

## `dgg_cells_shared_landscapes.csv.gz` — per-cell results

Gzip CSV, one row per domain cell (114,791,301), sorted by `seqnum`; 34 columns. Empty fields are
missing values. Pixel counts and `gaul0_code` are integers; `wrapped`, `popland` and `shland` are
1 (true) or 0 (false). Decimals written: 7 for coordinates and areas (about 1 cm and 0.1 m²),
5 for fractions (one 10 m pixel is ~0.000085 of a cell's land), 3 for distances, population and
densities. `popland` and `shland` were computed from the unrounded values, so a cell whose
fraction lies within 0.000005 of a threshold may appear to contradict its flag. Very small
values may be written in scientific notation (e.g. `8.5e-05`).

**Location and country**

| Column | Description |
|---|---|
| `seqnum` | cell id |
| `longitude`, `latitude` | cell centroid (degrees) |
| `gaul0_code`, `gaul0_name`, `iso3_code` | GAUL 2024 level 0 unit matched to the cell (empty if none) |
| `iso3_admin` | de facto country used in the country summaries (empty = unassigned) |
| `admin_source` | how `iso3_admin` was set: `gaul` (GAUL code as published), `crosswalk` (disputed GAUL feature, de facto administrator), `naturalearth` (Jammu & Kashmir / India–China border areas, per cell from Natural Earth), `naturalearth_remote` (> 50 km from GAUL, from Natural Earth), `unassigned` (Abyei, Bir Tawil, Spratly Islands), `beyond_max_dist` (> 50 km from GAUL and not in Natural Earth) |
| `continent` | continent from GAUL |
| `match_type` | `intersects` (centroid inside a GAUL polygon), `nearest` (nearest GAUL polygon), `ne_remote` (Natural Earth, remote islands) |
| `dist_km` | geodesic distance from the centroid to the matched polygon for `nearest` / `ne_remote` matches (km) |
| `cell_area_km2` | cell area (equal-area projection EPSG:6933) |
| `wrapped` | cell crosses the antimeridian (stored as east/west parts) |

**Population (WorldPop R2025A, 2024, 100 m constrained)**

| Column | Description |
|---|---|
| `pop_sum` | persons in the cell: sum of pixel values × fraction of each pixel covered by the cell (empty outside the raster's latitude extent) |
| `pop_density` | `pop_sum / cell_area_km2` |
| `valid_px_cov` | number of populated 100 m pixels covering the cell (coverage-weighted) |

**Earth Engine counts:** the `n_*` columns of `dgg_counts_2024.csv.gz` (above).

**Areas and classification**

| Column | Description |
|---|---|
| `land_area_km2` | land area: `n_land` × area of one 10 m pixel at the cell's latitude (WGS84 ellipsoid) |
| `habitat_area_km2` | habitat area, likewise from `n_habitat` |
| `nonhabitat_frac` | non-habitat fraction of the cell's land, `n_nonhabitat / n_land` (fraction, 0-1) |
| `habitat_frac` | habitat fraction of the cell's land, `n_habitat / n_land` (= `habitat_area_km2 / land_area_km2`; fraction, 0-1) |
| `pop_density_land` | persons per km² of land, `pop_sum / land_area_km2` (used for the classification) |
| `popland` | Populated Landscape: `n_land` > 0, `nonhabitat_frac` > 0.01 and `pop_density_land` ≥ 1 |
| `shland` | Shared Landscape: `popland` and `habitat_frac` ≥ 0.20 |

## Summary tables

`global_summary.csv` (one row, `group` = `Global`) and `country_summary.csv` (one row per de facto
country: `group` = name, `iso3`; `--` = unassigned) have these columns. Population fractions use
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
