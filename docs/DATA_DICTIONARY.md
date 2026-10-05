# Data dictionary

Grid: ISEA3H resolution 16 (DGGRID), WGS84 (EPSG:4326). Areas in km², population in persons.
Study domain of the counts and results: cells whose centroid latitude ≥ 60°S (114,791,301 cells).

## Grid

### `dgg_isea3h16_land_v4_shapefiles_partNNofMM.zip`

Parts of ~1 GB, each holding whole layers; together 474 zipped shapefiles, one per layer: `dgg_isea3h16_land_<n|s>LL_<e|w>LLL.zip` (10° × 10° tile
named by its south-west corner) and `dgg_isea3h16_land_south_cap.zip`. Each cell belongs to the
layer containing its centroid; together the layers hold every cell exactly once (126,509,398).

| Field | Type | Description |
|---|---|---|
| `seqnum` | text | cell id (DGGRID sequence number), stored as text because 9-digit integers overflow shapefile numeric fields |
| geometry | polygon / multipolygon | cell boundary; cells crossing ±180° are split into east and west parts |

### `dgg_isea3h16_land_v4_centroids_partNNofMM.zip`

Parts of ~1 GB, each holding whole files; together 474 CSV files, `<layer>_pts.csv`, one row per cell.

| Column | Type | Description |
|---|---|---|
| `seqnum` | integer | cell id |
| `longitude`, `latitude` | float | cell centroid, decimal degrees |

### `dgg_isea3h16_land_manifest.csv`

| Column | Description |
|---|---|
| `layer` | layer name |
| `lat0`, `lon0` | south-west corner of the tile (degrees) |
| `n_cells` | cells in the layer |
| `seconds` | build time of the layer |

## `dgg_counts_2024.parquet` — Earth Engine pixel counts

One row per domain cell. Counts are numbers of 10 m pixels (EPSG:4326, 10 m scale) whose centres
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

## `dgg_cells_shared_landscapes_partNNofMM.parquet` — per-cell results

One row per domain cell (114,791,301), split into parts of ~10 million consecutive rows ordered
by `seqnum` (read all parts together); 42 columns. Storage types: `seqnum`, `gaul0_code` int64;
text fields strings; `wrapped`, `popland`, `shland`, `has_gee`, `gee_geom_ok`
boolean; all other columns float64 (pixel counts are whole numbers stored as float64).

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
| `dist_deg`, `dist_km` | distance from the centroid to the matched polygon for `nearest` / `ne_remote` matches (degrees; geodesic km) |
| `cell_area_km2` | cell area (equal-area projection EPSG:6933) |
| `wrapped` | cell crosses the antimeridian (stored as east/west parts) |

**Population (WorldPop R2025A, 2024, 100 m constrained)**

| Column | Description |
|---|---|
| `pop_sum` | persons in the cell: sum of pixel values × fraction of each pixel covered by the cell (empty outside the raster's latitude extent) |
| `pop_density` | `pop_sum / cell_area_km2` |
| `valid_px_cov` | number of populated 100 m pixels covering the cell (coverage-weighted) |

**Earth Engine counts:** the `n_*` columns of `dgg_counts_2024.parquet` (above).

**Areas and classification**

| Column | Description |
|---|---|
| `land_area_km2` | land area: `n_land` × area of one 10 m pixel at the cell's latitude (WGS84 ellipsoid) |
| `habitat_area_km2` | habitat area, likewise from `n_habitat` |
| `land_pixel_km2`, `habitat_pixel_km2` | as above (pixel method; identical to `land_area_km2`, `habitat_area_km2`) |
| `land_baseline_km2`, `habitat_baseline_km2` | alternative cell-fraction areas, 1.18491 × `n_land` / `n_total` (resp. `n_habitat`); for comparison only |
| `nonhabitat_share` | non-habitat share of the cell's land, `n_nonhabitat / n_land` |
| `habitat_share` | habitat share of the cell's land, `n_habitat / n_land` (= `habitat_area_km2 / land_area_km2`) |
| `pop_density_land` | persons per km² of land, `pop_sum / land_area_km2` (used for the classification) |
| `popland` | Populated Landscape: `n_land` > 0, `nonhabitat_share` > 0.01 and `pop_density_land` ≥ 1 |
| `shland` | Shared Landscape: `popland` and `habitat_share` ≥ 0.20 |
| `has_gee` | Earth Engine counts present (true for every domain cell) |
| `gee_count_ratio` | QA: `n_total` × pixel area / `cell_area_km2` (≈ 1 for fully valid cells) |
| `gee_geom_ok` | QA: `gee_count_ratio` ≤ 1.02 (true for every cell) |

## Summary tables

`global_summary.csv`: `metric`, `area_km2`, `share` — total land; Populated Landscapes; total
habitat; habitat in Populated Landscapes; Shared Landscapes; Shared Landscapes within Populated
Landscapes. `share` is relative to total land, except for habitat in Populated Landscapes
(relative to total habitat) and Shared within Populated (relative to Populated Landscapes).

`country_summary.csv` (by de facto country: `group` = name, `iso3`; `--` = unassigned):

| Column | Description |
|---|---|
| `land_area_km2` | land |
| `populated_land_area_km2` | land in Populated Landscapes |
| `total_habitat_area_km2` | habitat |
| `habitat_in_populated_area_km2` | habitat in Populated Landscapes |
| `shl_land_area_km2` | land in Shared Landscapes |
| `shl_habitat_area_km2` | habitat in Shared Landscapes |
| `n_hex` | cells with land (`n_land` > 0) |
| `share_land_populated` | Populated Landscape share of land |
| `share_habitat_in_populated` | share of habitat in Populated Landscapes |
| `share_shl_land` | Shared Landscape share of land |
| `share_shl_within_populated` | Shared Landscape share of Populated Landscapes |
