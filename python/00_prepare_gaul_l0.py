"""
00_prepare_gaul_l0.py - GAUL 2024 level 0 (countries) from the FAO level 1 download.

FAO publishes GAUL 2024 at levels 1 and 2; level 0 is obtained by dissolving level 1 on the
country code, following the GEE community catalog tutorial (section 4, "Creating Level 0
Boundaries Using Python"):
https://gee-community-catalog.org/tutorials/examples/gaul_aggregate/#section-3-handling-gee-limitations

Input  (config.GAUL_L1_SHP): GAUL_2024_L1.shp from
       https://data.apps.fao.org/catalog/iso/34f97afc-6218-459a-971d-5af1162d318a
       (direct download: https://storage.googleapis.com/fao-maps-catalog-data/boundaries/GAUL_2024_L1.zip)
Output (config.GAUL_GPKG): gaul2024_l0.gpkg - one feature per gaul0_code; the other level 0
       attributes (gaul0_name, iso3_code, continent, ...) are taken from the first level 1 unit.

Run from the repository root:
    python -u python/00_prepare_gaul_l0.py
"""

import sys
from pathlib import Path

import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402


def main():
    src, dst = Path(config.GAUL_L1_SHP), Path(config.GAUL_GPKG)
    assert src.exists(), f"GAUL 2024 level 1 not found: {src} (see the docstring for the download)"
    gdf = gpd.read_file(src)
    gdf0 = gdf.dissolve(by="gaul0_code").drop(columns=["gaul1_code", "gaul1_name"])
    dst.parent.mkdir(parents=True, exist_ok=True)
    gdf0.to_file(dst, driver="GPKG")
    print(f"{src.name}: {len(gdf):,} level 1 units -> {len(gdf0):,} level 0 features -> {dst}")


if __name__ == "__main__":
    main()
