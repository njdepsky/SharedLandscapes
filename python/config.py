"""Configuration for the Python steps (04_download_worldpop.py, 05_postprocess.py).

All paths default to the repository's data/ folder; set the environment variable NRI_DATA_DIR
to use another location (e.g. a fast local disk), and NRI_WORK_DIR for the large intermediate
caches. This module must stay free of I/O: 05_postprocess.py's worker processes import it.
"""

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("NRI_DATA_DIR") or REPO / "data")
WORKDIR = Path(os.environ.get("NRI_WORK_DIR") or DATA_DIR / "work")   # caches (tens of GB)
OUT_DIR = DATA_DIR / "outputs"

# --- Inputs ---------------------------------------------------------------------------------
DGGDIR = DATA_DIR / "dgg"                                          # R/01_build_dgg.R output
GAUL_L1_SHP = DATA_DIR / "inputs" / "GAUL_2024_L1" / "GAUL_2024_L1.shp"   # FAO download (level 1)
GAUL_GPKG = DATA_DIR / "inputs" / "gaul2024_l0.gpkg"     # level 0, from 00_prepare_gaul_l0.py
NE_ADMIN0_SHP = DATA_DIR / "inputs" / "ne_10m_admin_0_countries" / "ne_10m_admin_0_countries.shp"

POP_YEAR = 2024
WORLDPOP_DIR = DATA_DIR / "inputs" / "worldpop" / f"R2025A_100m_{POP_YEAR}"
WORLDPOP = WORLDPOP_DIR / f"global_pop_{POP_YEAR}_CN_100m_R2025A_v1.vrt"   # 04_download_worldpop.py

# GEE counts (gee/03_extract_counts.js). Either ONE table (CSV or Parquet) - the default - or a
# FOLDER: then every file in it matching GEE_INPUT_PATTERN is read and combined.
GEE_INPUT = DATA_DIR / "gee_counts" / "dgg_counts_2024.parquet"
GEE_INPUT_PATTERN = "dgg_counts_*.csv"     # used only when GEE_INPUT is a folder

# --- Study domain ---------------------------------------------------------------------------
# Cells whose centroid latitude >= MIN_LAT (the DGG itself is global). Must match the layers
# extracted in Earth Engine (03_extract_counts.js: MIN_LAT).
MIN_LAT = -60.0

# --- Run control ----------------------------------------------------------------------------
N_WORKERS = 14              # parallel processes for the zonal population step
RERUN_STEP1 = False         # False: reuse the country-join cache (joins only new cells)
RERUN_STEP2 = False         # False: reuse the nearest-country cache when it matches
RERUN_ZONAL = False         # False: reuse per-layer population results
ID_FIELD = None             # cell-id field in the DGG shapefiles; None = auto-detect ("seqnum")

# --- Population -----------------------------------------------------------------------------
# WorldPop constrained products use nodata (never 0) over uninhabited land and water: count it
# as 0 people inside the raster's latitude extent (cells outside the extent get NaN).
NODATA_AS_ZERO = True

# --- Landscape classification (WL12-DWP40-WP100-SHL20) ---------------------------------------
LAND_AREA_METHOD = "pixel"   # "pixel": n_land x 10 m pixel area | "baseline": 1.18491 x n_land/n_total
NONHABITAT_THRESHOLD = 0.01            # Working Landscape: nonhabitat cover > 1% of land ...
POPULATION_DENSITY_THRESHOLD = 1.0     # ... and >= 1 person per km2 of land
SHL_HABITAT_THRESHOLD = 0.20           # Shared Landscape: Working + habitat cover >= 20%
POP_DENSITY_BASIS = "land"             # "land": people / land area | "cell": people / cell area
N_HEX_RULE = "n_land>0"                # cells counted in n_hex: "n_land>0" | "n_total>0" | "all"

# --- Country assignment ---------------------------------------------------------------------
GROUP_ISO3 = "iso3_admin"   # summaries by "iso3_admin" (de facto, below) or "iso3_code" (GAUL)
MAX_NEAREST_KM = 50.0       # unmatched cells farther than this from GAUL -> Natural Earth lookup
NE_SNAP_KM = 10.0           # ... accepted if inside a NE country or within this distance

# GAUL disputed features (lowercase "x" iso3 codes): de facto administrator, None = unassigned
ADMIN = {
    100: None,                                          # Abyei
    104: "FRA", 125: "FRA", 130: "FRA", 136: "FRA", 163: "FRA",   # Scattered Islands
    110: None,                                          # Bir Tawil
    133: "EGY",                                         # Hala'ib Triangle
    135: "KEN",                                         # Ilemi Triangle
    224: "CHN",                                         # Aksai Chin
    226: "IND",                                         # Arunachal Pradesh
    238: None,                                          # India-China border areas: per cell (NE)
    239: "ISR",                                         # Jerusalem-area 1949 no-man's-land strips
    245: None,                                          # Jammu and Kashmir: per cell (NE)
    249: "RUS",                                         # Kuril Islands
    262: "CHN",                                         # Paracel Islands
    267: "CHN",                                         # Scarborough Reef
    268: "JPN",                                         # Senkaku Islands
    270: None,                                          # Spratly Islands
    283: "GBR", 296: "GBR",                             # Akrotiri, Dhekelia
}
NE_CODES = [238, 245]                       # GAUL features resolved per cell with Natural Earth
NE_ALLOWED = ["IND", "PAK", "CHN", "KAS"]   # de facto administrators there ("KAS" = Siachen)
NE_RECODE = {"KAS": "IND"}
NE_ISO_RECODE = {                           # NE ADM0_A3 codes that differ from ISO3 / GAUL
    "SDS": "SSD", "PSX": "PSE", "SAH": "ESH",
    "CSI": "AUS",                           # Coral Sea Islands
    "BJN": "COL", "SER": "COL",             # Bajo Nuevo, Serranilla Bank (de facto)
}
