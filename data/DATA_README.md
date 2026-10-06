# data/ (not tracked by git)

Expected layout (override the root with the environment variable DATA_DIR):

    data/
      inputs/
        ne_10m_land/ne_10m_land.shp                      Natural Earth 10m land
        land-polygons-split-4326/land_polygons.shp       OSM land polygons (WGS84, split)
        ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp
        GAUL_2024_L1/GAUL_2024_L1.shp                    FAO GAUL 2024, level 1 (download)
        gaul2024_l0.gpkg                                 level 0, written by python/00_prepare_gaul_l0.py
        worldpop/R2025A_100m_2024/                       written by python/04_download_worldpop.py
      dgg/                                               written by R/01_build_dgg.R (or from Zenodo)
      gee_counts/                                        GEE CSVs, or dgg_counts_2024.parquet (Zenodo)
      work/                                              caches of python/05_postprocess.py (tens of GB)
      outputs/                                           per-cell tables and summaries/
