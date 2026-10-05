"""
05_postprocess.py - from GEE pixel counts to habitat and landscape extents and country summaries.

For every DGG cell in the study domain (centroid latitude >= MIN_LAT):
  1. country (GAUL 2024 level 0) by point-in-polygon of the cell centroid;
  2. nearest GAUL country for centroids outside every GAUL polygon (geodesic distance kept);
     beyond MAX_NEAREST_KM, Natural Earth 10m Admin 0 is used instead (remote islands);
  3. de facto country for disputed GAUL features (ADMIN crosswalk in config.py); Jammu &
     Kashmir and India-China border areas per cell from Natural Earth (IND / PAK / CHN);
  4. population: exact, coverage-weighted zonal sum of WorldPop (R2025A, 100 m, constrained)
     over each hexagon (exactextract), cached per DGG layer;
  5. per-cell table of country + population;
  6. join the GEE counts, classify, and write the per-cell result plus global and country
     summary tables.

Classification
  land_area_km2 (LAND_AREA_METHOD = "pixel") = n_land x area of one 10 m pixel (WGS84 ellipsoid,
      at the cell centroid latitude); habitat_area_km2 likewise with n_habitat.
      ("baseline": 1.18491 km2 x n_land / n_total - ocean pixels are not counted in n_total,
      so coastal cells are credited with nearly a full cell; kept for comparison.)
  habitat_share = habitat area / land area of the cell (= n_habitat / n_land).
  Populated Landscape (popland): n_land > 0, non-habitat share > 1% of land, and population
      density >= 1 person per km2 of land.
  Shared Landscape (shland): Populated Landscape with habitat share >= 20% of land.

Inputs (paths in config.py): DGG layers (R/01_build_dgg.R), GEE counts (gee/03_extract_counts.js;
one table, or every file matching a pattern in a folder), GAUL 2024 L0, Natural Earth 10m Admin 0,
WorldPop VRT (04_download_worldpop.py).

Outputs (config.OUT_DIR)
  dgg_cells_country_population.parquet   one row per domain cell: country fields, cell area,
                                         population (pop_sum, pop_density)
  dgg_cells_shared_landscapes.parquet    + GEE counts, land / habitat areas and shares, popland, shland
  summaries/global_summary.csv, country_summary.csv

Run from the repository root, in a terminal (worker processes cannot import a notebook):
    python -u python/05_postprocess.py 2>&1 | tee postprocess.log
"""

import os

os.environ.setdefault("GDAL_CACHEMAX", "1024")  # MB per process; set before GDAL initialises

import json
import multiprocessing as mp
import re
import shutil
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import dask.dataframe as dd
import dask_geopandas as dgpd
import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rasterio
from rasterio.windows import Window, from_bounds
import shapely
from exactextract import exact_extract
from shapely.affinity import translate
from pyproj import Geod
from shapely.geometry import box

# =============================================================================
# 0. CONFIG - user settings in config.py (imported by every worker process: no I/O here)
# =============================================================================
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import *  # noqa: E402,F401,F403  paths, domain, parameters, crosswalks

_TILE_LAT = re.compile(r"_([ns])(\d{2})_[ew]\d{3}$")


def layer_in_domain(stem):
    """DGG layers are owned by centroid latitude: a tile is in the domain iff lat0 >= MIN_LAT."""
    if stem.endswith("south_cap"):
        return MIN_LAT <= -90
    m = _TILE_LAT.search(stem)
    if not m:
        return True
    lat0 = int(m.group(2)) * (1 if m.group(1) == "n" else -1)
    return lat0 >= MIN_LAT


CRS = "EPSG:4326"
AREA_CRS = "EPSG:6933"  # WGS84 cylindrical equal-area, for measured cell areas
STEP1 = WORKDIR / "step1_country_join.parquet"            # incremental: base + increments
STEP1_INC = WORKDIR / "step1_country_join_increments"
STEP2 = WORKDIR / "step2_nearest_country.parquet"
ZONALDIR = WORKDIR / f"zonal_population_{Path(WORLDPOP).stem}"   # one parquet per DGG layer
OUT = OUT_DIR / "dgg_cells_country_population.parquet"
FINAL = OUT_DIR / "dgg_cells_shared_landscapes.parquet"
SUMMARY_DIR = OUT_DIR / "summaries"

HEX_AREA_KM2 = 1.18491  # nominal ISEA3H resolution-16 cell area (baseline land method)
GEE_PIXEL_DEG = 10 / 111319.49079327357  # GEE: 10 m scale in EPSG:4326 -> degrees per pixel
GEE_RATIO_MAX = 1.02  # QA: counted pixel area / cell area above this = invalid GEE geometry
WGS84_A_KM, WGS84_E2 = 6378.137, 0.00669437999014
GEE_ID = "seqnum16"
GEE_REQUIRED = ["n_total", "n_land", "n_nonhabitat", "n_habitat"]
GEOD = Geod(ellps="WGS84")
ID_CANDIDATES = ["seqnum", "SEQNUM", "global_id", "name", "id"]

# Reference cell area (QA only): ISEA3H is equal-area on the authalic sphere (DGGRID default);
# 10 x 3^r hexagon-equivalents cover it (the 12 pentagons are 5/6 of a hexagon each).
R_KM = 6371.007180918475
DGG_RES = 16
CELL_AREA_KM2 = 4 * np.pi * R_KM**2 / (10 * 3**DGG_RES)  # ~1.185 km2

COLS = ["gaul0_code", "gaul0_name", "iso3_code", "iso3_admin", "admin_source", "continent"]


# =============================================================================
# HELPERS (module level so spawned workers can import them)
# =============================================================================
def to_points(pdf):
    """pandas DataFrame with longitude/latitude -> GeoDataFrame of points."""
    return gpd.GeoDataFrame(
        pdf, geometry=gpd.points_from_xy(pdf["longitude"], pdf["latitude"]), crs=CRS
    )


def quiet_nearest(left, right, **kw):
    """sjoin_nearest in EPSG:4326 (distances in degrees; approximate by design)."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")
        return gpd.sjoin_nearest(left, right, **kw)


def _shift_east(coords):
    c = coords.copy()
    c[c[:, 0] < 0, 0] += 360.0
    return c


def fix_antimeridian(gs):
    """Split polygons whose vertices wrap across +/-180 (lon span > 180) into east/west parts."""
    b = gs.bounds
    wrapped = ((b["maxx"] - b["minx"]) > 180).to_numpy()
    if wrapped.any():
        fixed = []
        for g in gs[wrapped]:
            g2 = shapely.transform(g, _shift_east)
            east = g2.intersection(box(0, -90, 180, 90))
            west = translate(g2.intersection(box(180, -90, 360, 90)), xoff=-360)
            fixed.append(shapely.union(east, west))
        gs = gs.copy()
        gs.loc[wrapped] = gpd.GeoSeries(fixed, index=gs.index[wrapped], crs=gs.crs)
    return gs, wrapped


def _write_atomic(pdf, dst):
    """Write parquet via a hidden temp file + rename, so interrupted runs leave no partial file."""
    dst = Path(dst)
    tmp = dst.with_name(f".{dst.name}.tmp")
    pdf.to_parquet(tmp, index=False)
    os.replace(tmp, dst)


def zonal_tile(shp, out_dir, raster, id_field, rbounds):
    """Coverage-weighted WorldPop sum for every hexagon in one polygon tile -> parquet."""
    shp = Path(shp)
    dst = Path(out_dir) / f"{shp.stem}.parquet"
    g = gpd.read_file(shp, columns=[id_field])
    g = g.set_crs(CRS) if g.crs is None else g.to_crs(CRS)  # DGGRID SHPs may lack a .prj
    g["seqnum"] = pd.to_numeric(g[id_field]).astype("int64")
    g["geometry"], wrapped = fix_antimeridian(g.geometry)
    g["wrapped"] = wrapped
    g["cell_area_km2"] = g.geometry.to_crs(AREA_CRS).area / 1e6
    base = pd.DataFrame(g[["seqnum", "cell_area_km2", "wrapped"]])

    _, miny, _, maxy = g.total_bounds
    if maxy <= rbounds[1] or miny >= rbounds[3]:  # tile entirely outside raster lat extent
        _write_atomic(base.assign(pop_sum=0.0, valid_px_cov=0.0), dst)
        return shp.stem, len(base), int(wrapped.sum())

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        z = exact_extract(
            raster, g[["seqnum", "geometry"]], ["sum", "count"],
            include_cols=["seqnum"], output="pandas",
        )
    # "sum"   = sum(pixel value * coverage fraction) over valid (non-nodata) pixels
    # "count" = sum(coverage fraction) of valid pixels (pixel units)
    # robust to band-prefixed names (e.g. "band_1_sum") across exactextract versions
    ren = {c: ("pop_sum" if c.endswith("sum") else "valid_px_cov")
           for c in z.columns if c != "seqnum" and (c.endswith("sum") or c.endswith("count"))}
    assert len(ren) == 2, f"unexpected exactextract columns: {list(z.columns)}"
    z = z.rename(columns=ren).drop_duplicates("seqnum")
    _write_atomic(base.merge(z, on="seqnum", how="left"), dst)
    return shp.stem, len(base), int(wrapped.sum())


def run_zonal(tiles, id_field, rbounds, n_workers=N_WORKERS):
    """Run zonal_tile over tiles in parallel spawned processes, printing progress per tile."""
    t0, n, n_wrap, failed = time.time(), len(tiles), 0, []
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        futs = {
            ex.submit(zonal_tile, str(f), str(ZONALDIR), str(WORLDPOP), id_field, rbounds): f
            for f in tiles
        }
        for i, fut in enumerate(as_completed(futs), 1):
            f = futs[fut]
            try:
                stem, ncell, nw = fut.result()
                n_wrap += nw
                msg = f"{stem}: {ncell:,} cells"
            except Exception as e:  # keep going; failed tiles are simply not cached
                failed.append(f.name)
                msg = f"FAILED {f.name}: {type(e).__name__}: {e}"
            el = time.time() - t0
            print(f"  [{i}/{n}] {msg} | elapsed {el / 60:.1f} min | ETA {el / i * (n - i) / 60:.1f} min",
                  flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} tile(s) failed (re-run to retry): {failed}")
    return n_wrap


def subdivide(gdf, max_coords=2000, step=1.0):
    """Explode multipart features and cut large parts on a step-degree grid.

    The union of the pieces equals the original geometry, so point-to-feature distances are
    unchanged, but each distance computation now touches a small piece instead of a
    coastline with hundreds of thousands of vertices (the bottleneck of sjoin_nearest).
    """
    parts = gdf.explode(index_parts=False, ignore_index=True)
    ncoord = shapely.get_num_coordinates(np.asarray(parts.geometry))
    small, big = parts[ncoord <= max_coords], parts[ncoord > max_coords]
    if big.empty:
        return parts
    x0, y0, x1, y1 = big.total_bounds
    gx, gy = np.meshgrid(np.arange(np.floor(x0), np.ceil(x1), step),
                         np.arange(np.floor(y0), np.ceil(y1), step))
    gx, gy = gx.ravel(), gy.ravel()
    cells = shapely.box(gx, gy, gx + step, gy + step)
    geoms = np.asarray(big.geometry)
    ib, ic = shapely.STRtree(cells).query(geoms, predicate="intersects")
    # clip_by_rect only takes scalar bounds -> loop over grid cells, vectorised over the parts
    order = np.argsort(ic, kind="stable")
    ib, ic = ib[order], ic[order]
    ucells, starts = np.unique(ic, return_index=True)
    ends = np.append(starts[1:], len(ic))
    pieces = np.empty(len(ib), dtype=object)
    for c, s0, s1 in zip(ucells, starts, ends):
        x, y = float(gx[c]), float(gy[c])
        pieces[s0:s1] = shapely.clip_by_rect(geoms[ib[s0:s1]], x, y, x + step, y + step)
    keep = ~shapely.is_empty(pieces)
    attrs = big.drop(columns=big.geometry.name).iloc[ib[keep]].reset_index(drop=True)
    cut = gpd.GeoDataFrame(attrs, geometry=pieces[keep], crs=gdf.crs)
    return pd.concat([small, cut], ignore_index=True)


def geodesic_km(points, geoms):
    """WGS84 geodesic distance (km) from each point to the closest point of its geometry.

    The closest point is found in lon/lat space (shapely), then measured geodesically -
    accurate enough for a tens-of-km threshold; 0 for points inside their polygon.
    """
    lines = shapely.shortest_line(np.asarray(points), np.asarray(geoms))
    c = shapely.get_coordinates(lines).reshape(-1, 2, 2)
    _, _, d = GEOD.inv(c[:, 0, 0], c[:, 0, 1], c[:, 1, 0], c[:, 1, 1])
    return d / 1000.0


def raster_total(path, rows=1024):
    """Sum of all valid pixels in the raster (QA only), cached per raster file.

    Reads full-width strips of `rows` rows (a 100 m VRT has millions of small blocks, so
    block-by-block reading takes hours). The result is cached in WORKDIR/raster_totals.json,
    keyed by path + size + modification time, so it is computed once per raster.
    """
    p = Path(path).resolve()
    st = p.stat()
    key = f"{p}|{st.st_size}|{int(st.st_mtime)}"
    cache = WORKDIR / "raster_totals.json"
    data = json.loads(cache.read_text()) if cache.exists() else {}
    if key in data:
        return data[key]
    t0 = time.time()
    tot = 0.0
    with rasterio.open(p) as src:
        nod = src.nodata
        for r0 in range(0, src.height, rows):
            a = src.read(1, window=Window(0, r0, src.width, min(rows, src.height - r0)))
            m = np.isfinite(a) if nod is None else ((a != nod) & np.isfinite(a))
            tot += float(a[m].sum(dtype="float64"))
    data[key] = tot
    cache.write_text(json.dumps(data, indent=1))
    print(f"(raster total computed in {(time.time() - t0) / 60:.1f} min; cached for later runs)", flush=True)
    return tot


# =============================================================================
# GEE COUNTS: one table, or every file matching GEE_INPUT_PATTERN in a folder
# =============================================================================
def gee_input_files():
    p = Path(GEE_INPUT)
    if p.is_dir():
        files = sorted(f for f in p.glob(GEE_INPUT_PATTERN) if f.is_file())
        assert files, f"no files matching {GEE_INPUT_PATTERN!r} in {p}"
        return files
    assert p.exists(), f"GEE input not found: {p}"
    return [p]


def gee_columns(f):
    f = Path(f)
    return pq.read_schema(f).names if f.suffix == ".parquet" else pd.read_csv(f, nrows=0).columns.tolist()


def load_gee_counts():
    """All GEE count files -> one frame (seqnum + n_* bands). A cell may appear in several files
    only with identical values (it is then kept once); conflicting repeats stop the run."""
    files = gee_input_files()
    n_cols = [c for c in gee_columns(files[0]) if c.startswith("n_")]
    parts = []
    for f in files:
        d = (pd.read_parquet(f, columns=[GEE_ID] + n_cols) if f.suffix == ".parquet"
             else pd.read_csv(f, usecols=[GEE_ID] + n_cols, dtype={GEE_ID: "string"}))
        x = pd.to_numeric(d[GEE_ID]).to_numpy(dtype="float64")
        assert np.isfinite(x).all() and (x == np.floor(x)).all(), f"{f.name}: non-integer {GEE_ID}"
        d[GEE_ID] = x.astype("int64")
        for c in n_cols:
            d[c] = pd.to_numeric(d[c]).astype("float64")
        parts.append(d)
    gee = pd.concat(parts, ignore_index=True).rename(columns={GEE_ID: "seqnum"})
    dup = gee["seqnum"].duplicated(keep=False)
    if dup.any():
        n_conf = int((gee.loc[dup].groupby("seqnum")[n_cols].nunique().max(axis=1) > 1).sum())
        assert n_conf == 0, f"GEE input: {n_conf:,} cells repeated with different counts"
        gee = gee.drop_duplicates("seqnum")
    print(f"GEE counts: {len(files)} file(s), {len(gee):,} cells, bands {n_cols}", flush=True)
    return gee


# =============================================================================
# STEP 6 - GEE JOIN, LANDSCAPE CLASSIFICATION, SUMMARY TABLES
# =============================================================================
SUMMARY_COLS = [
    "land_area_km2", "populated_land_area_km2",
    "total_habitat_area_km2", "habitat_in_populated_area_km2",
    "shl_land_area_km2", "shl_habitat_area_km2", "n_hex",
    "share_land_populated", "share_habitat_in_populated",
    "share_shl_land", "share_shl_within_populated",
]
OUT_STR_COLS = ["gaul0_name", "iso3_code", "iso3_admin", "admin_source", "continent", "match_type"]


def _num(s):
    """Series (numpy-, nullable- or arrow-backed) -> float64 numpy array, NaN for missing."""
    return s.to_numpy(dtype="float64", na_value=np.nan)


def _div(num, den):
    """Element-wise num/den; NaN where den is 0 or not finite."""
    num, den = np.asarray(num, dtype="float64"), np.asarray(den, dtype="float64")
    out = np.full(np.broadcast(num, den).shape, np.nan)
    np.divide(num, den, out=out, where=np.isfinite(den) & (den != 0))
    return out


def pixel_area_km2(lat_deg):
    """WGS84 ellipsoidal area (km2) of one GEE_PIXEL_DEG x GEE_PIXEL_DEG pixel at latitude."""
    phi = np.radians(lat_deg)
    s2 = np.sin(phi) ** 2
    m = WGS84_A_KM * (1 - WGS84_E2) / (1 - WGS84_E2 * s2) ** 1.5  # meridional radius
    n = WGS84_A_KM / np.sqrt(1 - WGS84_E2 * s2)                    # prime-vertical radius
    return np.radians(GEE_PIXEL_DEG) ** 2 * m * n * np.cos(phi)


def classify(df):
    """Add landscape fields (see module docstring) to a frame with GEE + pop columns."""
    n_total, n_land, n_non, n_hab = (_num(df[c]) for c in GEE_REQUIRED)
    px = pixel_area_km2(_num(df["latitude"]))
    # GEE geometry check: counted pixel area / real cell area. ~1 for sound cells; >> 1 where GEE
    # read an antimeridian-crossing hexagon as a globe-spanning band (counts are then invalid).
    ratio = _div(n_total * px, _num(df["cell_area_km2"]))
    with np.errstate(invalid="ignore"):
        geom_ok = ~(ratio > GEE_RATIO_MAX)
    df["gee_count_ratio"] = ratio
    df["gee_geom_ok"] = geom_ok
    base_land, base_hab = _div(n_land, n_total) * HEX_AREA_KM2, _div(n_hab, n_total) * HEX_AREA_KM2
    variants = {
        "baseline": (base_land, base_hab),
        # pixel area where GEE geometry is sound; baseline fallback on flagged cells (see QA)
        "pixel": (np.where(geom_ok, n_land * px, base_land), np.where(geom_ok, n_hab * px, base_hab)),
    }
    if LAND_AREA_METHOD not in variants:
        raise ValueError(f"LAND_AREA_METHOD must be 'baseline' or 'pixel', not {LAND_AREA_METHOD!r}")
    for k, (lv, hv) in variants.items():
        df[f"land_{k}_km2"] = lv
        df[f"habitat_{k}_km2"] = hv
    land, habitat = variants[LAND_AREA_METHOD]
    nonhab_cover = _div(n_non, n_land)
    hab_cover = _div(n_hab, n_land)
    pop = _num(df["pop_sum"])
    if POP_DENSITY_BASIS == "land":
        dens = _div(pop, land)
    elif POP_DENSITY_BASIS == "cell":
        dens = _div(pop, _num(df["cell_area_km2"]))
    else:
        raise ValueError(f"POP_DENSITY_BASIS must be 'land' or 'cell', not {POP_DENSITY_BASIS!r}")
    with np.errstate(invalid="ignore"):  # NaN comparisons -> False (not populated)
        popland = (n_land > 0) & (nonhab_cover > NONHABITAT_THRESHOLD) & (
            dens >= POPULATION_DENSITY_THRESHOLD)
        shland = popland & (hab_cover >= SHL_HABITAT_THRESHOLD)
    df["land_area_km2"] = land
    df["habitat_area_km2"] = habitat
    df["nonhabitat_share"] = nonhab_cover   # non-habitat share of the cell's land
    df["habitat_share"] = hab_cover         # habitat share of the cell's land (= habitat / land area)
    df["pop_density_land"] = dens
    df["popland"] = popland
    df["shland"] = shland
    df["has_gee"] = np.isfinite(n_total)
    return df


def _finish(g):
    """Aggregated sums (land, habitat, p_land, p_hab, s_land, s_hab, n_hex) -> summary columns."""
    r = pd.DataFrame(index=g.index)
    r["land_area_km2"] = g["land"]
    r["populated_land_area_km2"] = g["p_land"]
    r["total_habitat_area_km2"] = g["habitat"]
    r["habitat_in_populated_area_km2"] = g["p_hab"]
    r["shl_land_area_km2"] = g["s_land"]
    r["shl_habitat_area_km2"] = g["s_hab"]
    r["n_hex"] = g["n_hex"].astype("int64")
    r["share_land_populated"] = _div(r["populated_land_area_km2"], r["land_area_km2"])
    r["share_habitat_in_populated"] = _div(r["habitat_in_populated_area_km2"], r["total_habitat_area_km2"])
    r["share_shl_land"] = _div(r["shl_land_area_km2"], r["land_area_km2"])
    r["share_shl_within_populated"] = _div(r["shl_land_area_km2"], r["populated_land_area_km2"])
    return r


def _write_csv(frame, path):
    tmp = Path(path).with_name(f".{Path(path).name}.tmp")
    frame.to_csv(tmp, index=False, float_format="%.4f")
    os.replace(tmp, path)


def summarize(final, iso3_to_name):
    """Global and country summary tables over cells with GEE data."""
    has = final["has_gee"].to_numpy(bool)
    land = np.nan_to_num(_num(final["land_area_km2"]))[has]
    hab = np.nan_to_num(_num(final["habitat_area_km2"]))[has]
    p = final["popland"].to_numpy(bool)[has]
    s = final["shland"].to_numpy(bool)[has]
    if N_HEX_RULE == "all":
        cnt = np.ones(int(has.sum()), dtype="int64")
    elif N_HEX_RULE in ("n_land>0", "n_total>0"):
        col = "n_land" if N_HEX_RULE == "n_land>0" else "n_total"
        cnt = (np.nan_to_num(_num(final[col]))[has] > 0).astype("int64")
    else:
        raise ValueError(f"N_HEX_RULE must be 'n_land>0', 'n_total>0' or 'all', not {N_HEX_RULE!r}")
    base = pd.DataFrame({"land": land, "habitat": hab, "p_land": land * p, "p_hab": hab * p,
                         "s_land": land * s, "s_hab": hab * s, "cnt": cnt})
    aggs = {c: (c, "sum") for c in ["land", "habitat", "p_land", "p_hab", "s_land", "s_hab"]}
    aggs["n_hex"] = ("cnt", "sum")

    def by(col):
        c = final[col]
        c = c if isinstance(c.dtype, pd.CategoricalDtype) else c.astype("category")
        base["_k"] = pd.Categorical.from_codes(c.cat.codes.to_numpy()[has], c.cat.categories)
        g = base.groupby("_k", observed=True, dropna=False).agg(**aggs).reset_index()
        base.drop(columns="_k", inplace=True)
        return g

    # --- global
    T = base[["land", "habitat", "p_land", "p_hab", "s_land", "s_hab", "cnt"]].sum()
    glob = pd.DataFrame([
        ("Total land area", T.land, 1.0),
        ("Populated Landscapes", T.p_land, _div(T.p_land, T.land).item()),
        ("Total habitat area", T.habitat, _div(T.habitat, T.land).item()),
        ("Habitat in Populated Landscapes", T.p_hab, _div(T.p_hab, T.habitat).item()),
        ("Shared Landscapes (SHL)", T.s_land, _div(T.s_land, T.land).item()),
        ("SHL within Populated Landscapes", T.s_land, _div(T.s_land, T.p_land).item()),
    ], columns=["metric", "area_km2", "share"])

    # --- country
    g = by(GROUP_ISO3)
    names = {str(k).upper(): v for k, v in iso3_to_name.items()}
    keys = [None if pd.isna(k) else str(k).upper() for k in g["_k"]]
    country = _finish(g)
    country.insert(0, "group", ["Unassigned" if k is None else names.get(k, k) for k in keys])
    country.insert(1, "iso3", ["--" if k is None else k for k in keys])
    # deterministic order: land area descending, ties broken by iso3 (stable sort)
    country = country.sort_values(["land_area_km2", "iso3"], ascending=[False, True],
                                  kind="mergesort")[["group", "iso3"] + SUMMARY_COLS]
    return glob, country, int(T.cnt)


def run_step6(iso3_to_name):
    t0 = time.time()
    print("\n" + "=" * 78 + "\nStep 6: GEE join, landscape classification, summaries", flush=True)

    # --- DGG cell table from step 5 (strings dictionary-encoded -> pandas categoricals)
    cells = pq.read_table(OUT, read_dictionary=OUT_STR_COLS).to_pandas()
    cells["seqnum"] = cells["seqnum"].to_numpy(dtype="int64")
    # Fixed cell order (by id): the per-cell output and every floating-point sum below are then
    # identical however the step 1-4 caches were built (from scratch or incrementally).
    cells = cells.sort_values("seqnum", kind="mergesort", ignore_index=True)

    # --- GEE counts (all n_* bands carried through)
    gee = load_gee_counts()
    in_cells = gee["seqnum"].isin(cells["seqnum"])
    print(f"Step 6: GEE rows {len(gee):,} | not in DGG cell table: {(~in_cells).sum():,} "
          f"(expect 0; dropped if any)", flush=True)
    gee = gee[in_cells]

    final = cells.merge(gee, on="seqnum", how="left")
    del cells, gee
    final = classify(final)
    has = final["has_gee"].to_numpy(bool)
    lat = _num(final["latitude"])
    print(f"Step 6: DGG cells {len(final):,} | with GEE data {has.sum():,} | without {(~has).sum():,} "
          f"(of which lat >= -60: {((~has) & (lat >= -60)).sum():,}; expect 0)", flush=True)
    nt, nl = _num(final["n_total"]), _num(final["n_land"])
    print(f"Step 6: n_total == 0: {(has & (nt == 0)).sum():,} cells (contribute 0 land) | "
          f"n_land > n_total: {(has & (nl > nt)).sum():,} (expect 0)", flush=True)
    print(f"Step 6: Populated Landscape cells {final['popland'].sum():,} | Shared Landscape cells "
          f"{final['shland'].sum():,} "
          f"(density basis: {POP_DENSITY_BASIS}; land area method: {LAND_AREA_METHOD})", flush=True)
    print(f"Step 6: total land - baseline {np.nansum(_num(final['land_baseline_km2'])):,.1f} km2 | "
          f"pixel {np.nansum(_num(final['land_pixel_km2'])):,.1f} km2", flush=True)
    bad = has & ~final["gee_geom_ok"].to_numpy(bool)
    print(f"Step 6: GEE geometry flags (count ratio > {GEE_RATIO_MAX}): {bad.sum():,} cells "
          f"({final['wrapped'].to_numpy(bool)[bad].sum():,} antimeridian-wrapped) | baseline land in them "
          f"{np.nansum(_num(final['land_baseline_km2'])[bad]):,.1f} km2 (their counts are invalid; "
          f"pixel method falls back to baseline there)", flush=True)

    tmp = FINAL.with_name(f".{FINAL.name}.tmp")
    final.to_parquet(tmp, index=False, row_group_size=2_000_000)
    os.replace(tmp, FINAL)
    print(f"Step 6: wrote {FINAL} ({len(final):,} rows)", flush=True)

    glob, country, n_hex_total = summarize(final, iso3_to_name)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    stale = SUMMARY_DIR / "region_summary.csv"   # no longer produced
    if stale.exists():
        stale.unlink()
    _write_csv(glob, SUMMARY_DIR / "global_summary.csv")
    _write_csv(country, SUMMARY_DIR / "country_summary.csv")

    # --- QA: groups must close on the global totals
    tl = float(glob.loc[glob.metric == "Total land area", "area_km2"].iloc[0])
    print(f"\nStep 6 QA: land total {tl:,.4f} km2 | country sum {country.land_area_km2.sum():,.4f} | "
          f"n_hex ({N_HEX_RULE}) global {n_hex_total:,} | country {int(country.n_hex.sum()):,}")
    ua = country[country.iso3 == "--"]
    if len(ua):
        print(f"Step 6 QA: unassigned country cells {int(ua.n_hex.iloc[0]):,} "
              f"({ua.land_area_km2.iloc[0]:,.1f} km2 land)")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print("\nglobal_summary.csv:\n" + glob.to_string(index=False))
        print(f"\ncountry_summary.csv ({len(country)} rows), top 10:\n"
              + country.iloc[:10, :8].to_string(index=False))
    print(f"\nStep 6: done ({(time.time() - t0) / 60:.1f} min) -> {SUMMARY_DIR}", flush=True)


# =============================================================================
# MAIN PIPELINE
# =============================================================================
def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ZONALDIR.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. INPUTS
    # -------------------------------------------------------------------------
    # --- GAUL 2024 Level 0
    gaul_l0 = gpd.read_file(GAUL_GPKG)

    # --- Natural Earth 10m Admin 0 Countries (de facto boundaries; Kashmir region only).
    #     ADM0_A3 rather than ISO_A3: NE's ISO_A3 is "-99" for some countries.
    ne = gpd.read_file(NE_ADMIN0_SHP)[["ADM0_A3", "geometry"]].to_crs(CRS)

    # --- WorldPop (people per pixel): sanity checks
    with rasterio.open(WORLDPOP) as src:
        assert src.crs is not None and src.crs.to_epsg() == 4326, f"WorldPop CRS is {src.crs}"
        print(f"WorldPop file: {WORLDPOP.name}")
        print(f"WorldPop: {src.width}x{src.height}, res={src.res}, bounds={src.bounds}, "
              f"dtype={src.dtypes[0]}, nodata={src.nodata}")
        print("WorldPop tags:", src.tags())
        rbounds = tuple(src.bounds)  # (left, bottom, right, top)
        # pre-flight guard: a raster without declared nodata (e.g. a mis-built VRT mosaic) would
        # feed fill values into the zonal sums; a populated control window must look populated
        assert src.nodata is not None, "WorldPop raster declares no nodata value - fix/rebuild it"
        ctl = src.read(1, window=from_bounds(80, 25, 84, 27, src.transform))  # Gangetic plain
        ctl_valid = ctl[ctl != src.nodata] if not np.isnan(src.nodata) else ctl[~np.isnan(ctl)]
        print(f"WorldPop control window (Gangetic plain): populated {np.mean(ctl > 0):.1%}, "
              f"negative valid values {(ctl_valid < 0).sum():,}")
        assert np.mean(ctl > 0) > 0.3 and not (ctl_valid < 0).any(), "WorldPop control window check failed"

    # --- DGG polygon tiles (10x10 deg SHPs)
    shp_all = sorted(DGGDIR.glob("*.shp"))
    shp_tiles = [f for f in shp_all if layer_in_domain(f.stem)]
    print(f"Domain: centroid latitude >= {MIN_LAT:g} -> {len(shp_tiles)} of {len(shp_all)} layers")
    assert shp_tiles, f"no .shp tiles found in {DGGDIR}"
    t = gpd.read_file(shp_tiles[0], rows=5)
    id_field = ID_FIELD or next((c for c in ID_CANDIDATES if c in t.columns), None)
    assert id_field, f"no cell-ID field found; columns = {list(t.columns)} -> set ID_FIELD"
    print(f"DGG polygon tiles: {len(shp_tiles)} | ID field: '{id_field}' | "
          f"columns: {list(t.columns)} | CRS: {t.crs}")

    # --- GEE counts (used in step 6): fail early if absent or incomplete
    gee_files = gee_input_files()
    gee_cols = set(gee_columns(gee_files[0]))
    missing_gee = [c for c in [GEE_ID] + GEE_REQUIRED if c not in gee_cols]
    assert not missing_gee, f"GEE input lacks columns {missing_gee}; has {sorted(gee_cols)}"
    print(f"GEE input: {len(gee_files)} file(s) ({Path(GEE_INPUT).name})")

    # --- DGG cell centroids (one CSV per 10-degree tile).
    #     seqnum is read as float64 (robust to ids written in scientific notation by some
    #     tools); the cast back to int64 is exact (no NaN, integral, << 2**53).
    pts_files = [str(f) for f in sorted(DGGDIR.glob("*_pts.csv"))
                 if layer_in_domain(f.name[: -len("_pts.csv")])]
    dgg = dd.read_csv(
        pts_files,
        dtype={"seqnum": "float64", "longitude": "float64", "latitude": "float64"},
    )
    dgg["seqnum"] = dgg["seqnum"].astype("int64")

    dgg_gdf = dgpd.from_dask_dataframe(
        dgg,
        geometry=dgpd.points_from_xy(dgg, x="longitude", y="latitude", crs=CRS),
    ).set_crs(CRS, allow_override=True)

    # -------------------------------------------------------------------------
    # 2. GAUL PREP + DISPUTED-AREA CROSSWALK
    # -------------------------------------------------------------------------
    gaul = gaul_l0.to_crs(CRS).reset_index(drop=True)
    assert gaul["gaul0_code"].is_unique, "gaul0_code not unique - attribute join would duplicate rows"

    isx = gaul["iso3_code"].str.startswith("x")
    missing = set(gaul.loc[isx, "gaul0_code"]) - set(ADMIN)
    assert not missing, f"x-code features missing from ADMIN: {missing}"

    gaul["iso3_admin"] = np.where(isx, gaul["gaul0_code"].map(ADMIN), gaul["iso3_code"])
    gaul["admin_source"] = np.select(
        [~isx, gaul["iso3_admin"].notna()], ["gaul", "crosswalk"], default="unassigned"
    )
    gaul_key = gaul[["gaul0_code", "geometry"]]  # geometry + key only, for spatial joins
    gaul_attr = pd.DataFrame(gaul[COLS])  # attributes, joined on gaul0_code later
    gaul_attr["gaul0_code"] = gaul_attr["gaul0_code"].astype("Int64")

    # ISO3 -> GAUL feature, for mapping Natural Earth countries back to GAUL (non-disputed only).
    # Where one ISO3 has several GAUL features (e.g. AUS: Australia + Ashmore and Cartier Is.;
    # PYF: French Polynesia + Clipperton I.), use the LARGEST feature (the main country).
    g_nx = gaul.loc[~isx, ["iso3_code", "gaul0_code", "gaul0_name", "geometry"]].copy()
    g_nx["area"] = g_nx.geometry.to_crs(AREA_CRS).area
    g_nx = g_nx.sort_values("area", ascending=False)
    dup = g_nx[g_nx["iso3_code"].duplicated(keep=False)]
    if len(dup):
        chosen = dup.drop_duplicates("iso3_code")[["iso3_code", "gaul0_name"]].values.tolist()
        print(f"note: ISO3 codes with >1 GAUL feature; largest used for NE mapping: {chosen}")
    iso3_to_gaul = g_nx.drop_duplicates("iso3_code").set_index("iso3_code")["gaul0_code"]

    # ISO3 (incl. x-codes) -> display name for country summaries: largest GAUL feature per code
    g_all = gaul[["iso3_code", "gaul0_name", "geometry"]].copy()
    g_all["area"] = g_all.geometry.to_crs(AREA_CRS).area
    iso3_to_name = (g_all.sort_values("area", ascending=False)
                    .drop_duplicates("iso3_code").set_index("iso3_code")["gaul0_name"])

    # -------------------------------------------------------------------------
    # 3. STEP 1 - POINT-IN-POLYGON (expensive; run once and persist)
    # -------------------------------------------------------------------------
    def assign_intersects(part):
        part = part.set_crs(CRS, allow_override=True)  # guard against CRS dropped per partition
        j = gpd.sjoin(part, gaul_key, how="left", predicate="intersects")
        # points on shared borders / in overlapping polygons match twice -> keep one
        j = j.drop_duplicates(subset="seqnum", keep="first").drop(columns="index_right")
        j["gaul0_code"] = j["gaul0_code"].astype("Int64")  # consistent dtype across partitions
        return j

    def step1_files():
        return sorted(STEP1.glob("*.parquet")) + sorted(STEP1_INC.glob("inc_*.parquet/*.parquet"))

    meta1 = assign_intersects(dgg_gdf._meta)
    if RERUN_STEP1 or not step1_files():
        print("Step 1: running spatial join for all cells ...", flush=True)
        if STEP1_INC.exists():
            shutil.rmtree(STEP1_INC)
        joined = dgg_gdf.map_partitions(assign_intersects, meta=meta1)
        joined.to_parquet(STEP1, write_index=False, overwrite=True)
        stale = np.array([], dtype="int64")
    else:
        grid_ids = np.unique(dgg["seqnum"].compute().to_numpy())
        cached = np.unique(dd.read_parquet(step1_files(), columns=["seqnum"])["seqnum"]
                           .compute().to_numpy().astype("int64"))
        missing = np.setdiff1d(grid_ids, cached, assume_unique=True)
        stale = np.setdiff1d(cached, grid_ids, assume_unique=True)
        print(f"Step 1: grid {len(grid_ids):,} cells | cached {len(cached):,} | to join {len(missing):,} "
              f"| cached but not in grid (ignored) {len(stale):,}", flush=True)
        if len(missing):
            miss_idx = pd.Index(missing)
            sub = dgg_gdf.map_partitions(lambda p: p[p["seqnum"].isin(miss_idx)], meta=dgg_gdf._meta)
            STEP1_INC.mkdir(parents=True, exist_ok=True)
            inc = STEP1_INC / f"inc_{len(list(STEP1_INC.glob('inc_*.parquet'))) + 1:03d}.parquet"
            t0 = time.time()
            sub.map_partitions(assign_intersects, meta=meta1).to_parquet(inc, write_index=False)
            print(f"Step 1: joined {len(missing):,} new cells -> {inc.name} "
                  f"({(time.time() - t0) / 60:.1f} min)", flush=True)

    # Plain (non-geo) read: geometry is rebuilt from lon/lat only where needed
    j1 = dd.read_parquet(step1_files(), columns=["seqnum", "longitude", "latitude", "gaul0_code"])
    j1["gaul0_code"] = j1["gaul0_code"].astype("Int64")
    if len(stale):  # cells from an earlier grid version that this grid no longer contains
        stale_idx = pd.Index(stale)
        j1 = j1.map_partitions(lambda p: p[~p["seqnum"].isin(stale_idx)], meta=j1._meta)

    # -------------------------------------------------------------------------
    # 4. STEP 2 - NEAREST GAUL FEATURE FOR UNMATCHED CELLS
    # -------------------------------------------------------------------------
    um = j1[j1["gaul0_code"].isna()][["seqnum", "longitude", "latitude"]].compute()
    if STEP2.exists() and not RERUN_STEP2:
        cached2 = pd.read_parquet(STEP2, columns=["seqnum"])["seqnum"]
        if set(cached2.tolist()) != set(um["seqnum"].tolist()):
            print(f"Step 2: cached nearest set ({len(cached2):,}) differs from unmatched cells "
                  f"({len(um):,}) - recomputing", flush=True)
            STEP2.unlink()
    if RERUN_STEP2 or not STEP2.exists():
        t0 = time.time()
        print(f"Step 2: unmatched cells = {len(um):,}", flush=True)

        gaul_near = subdivide(gaul_key)  # same union, small pieces -> fast nearest queries
        print(f"Step 2: GAUL subdivided into {len(gaul_near):,} pieces "
              f"({time.time() - t0:.0f} s)", flush=True)

        nn = quiet_nearest(to_points(um), gaul_near, how="left", distance_col="dist_deg")
        nn = nn.sort_values(["seqnum", "dist_deg"]).drop_duplicates("seqnum")  # closest piece
        nn["dist_km"] = geodesic_km(
            nn.geometry.values, gaul_near.geometry.values[nn["index_right"].to_numpy()]
        )
        near = pd.DataFrame(nn.drop(columns=["geometry", "index_right"])).assign(match_type="nearest")
        near["gaul0_code"] = near["gaul0_code"].astype("Int64")

        # --- Natural Earth rescue for cells farther than MAX_NEAREST_KM from any GAUL feature
        far = near["dist_km"] > MAX_NEAREST_KM
        print(f"Step 2: {far.sum():,} cells > {MAX_NEAREST_KM:g} km from GAUL -> Natural Earth lookup",
              flush=True)
        if far.any():
            ne_sub = subdivide(ne)
            nn2 = quiet_nearest(to_points(near.loc[far, ["seqnum", "longitude", "latitude"]]),
                                ne_sub, how="left", distance_col="d_deg")
            nn2 = nn2.sort_values(["seqnum", "d_deg"]).drop_duplicates("seqnum")
            nn2["ne_km"] = geodesic_km(
                nn2.geometry.values, ne_sub.geometry.values[nn2["index_right"].to_numpy()]
            )
            code = nn2["ADM0_A3"].replace(NE_ISO_RECODE).map(iso3_to_gaul)
            ok = (nn2["ne_km"] <= NE_SNAP_KM) & code.notna()
            upd = pd.DataFrame({"seqnum": nn2.loc[ok, "seqnum"].to_numpy(),
                                "gaul0_code_ne": code[ok].astype("Int64").to_numpy(),
                                "dist_km_ne": nn2.loc[ok, "ne_km"].to_numpy()})
            near = near.merge(upd, on="seqnum", how="left")
            hit = near["gaul0_code_ne"].notna()
            near.loc[hit, "gaul0_code"] = near.loc[hit, "gaul0_code_ne"]
            near.loc[hit, "dist_km"] = near.loc[hit, "dist_km_ne"]
            near.loc[hit, "dist_deg"] = np.nan
            near.loc[hit, "match_type"] = "ne_remote"
            near = near.drop(columns=["gaul0_code_ne", "dist_km_ne"])
            unmapped = nn2.loc[(nn2["ne_km"] <= NE_SNAP_KM) & code.isna(), "ADM0_A3"].value_counts()
            print(f"Step 2: NE rescued {hit.sum():,} cells; still unresolved "
                  f"{far.sum() - hit.sum():,}", flush=True)
            if len(unmapped):
                print("Step 2: NE countries hit but not mappable to GAUL ISO3 "
                      f"(left unassigned): {unmapped.to_dict()}", flush=True)

        _write_atomic(near, STEP2)
        print(f"Step 2: nearest join done ({(time.time() - t0) / 60:.1f} min)", flush=True)
    else:
        near = pd.read_parquet(STEP2)
        print(f"Step 2: reusing {STEP2} ({len(near):,} cells)", flush=True)

    # -------------------------------------------------------------------------
    # 5. STEP 3 - NATURAL EARTH PER-CELL DE FACTO ASSIGNMENT (Kashmir region)
    # -------------------------------------------------------------------------
    kas_pts = pd.concat(
        [
            j1[j1["gaul0_code"].isin(NE_CODES)][["seqnum", "longitude", "latitude"]].compute(),
            near.loc[near["gaul0_code"].isin(NE_CODES), ["seqnum", "longitude", "latitude"]],
        ],
        ignore_index=True,
    )
    kas = to_points(kas_pts)
    ne_k = ne[ne["ADM0_A3"].isin(NE_ALLOWED)]

    k = gpd.sjoin(kas, ne_k, how="left", predicate="intersects").drop_duplicates("seqnum")
    miss = k["ADM0_A3"].isna()
    if miss.any():  # NE gaps/edges or outside allowed countries -> nearest allowed country
        k_nn = quiet_nearest(kas[kas["seqnum"].isin(k.loc[miss, "seqnum"])], ne_k, how="left")
        k = pd.concat([k[~miss], k_nn.drop_duplicates("seqnum")])

    kas_admin = k.set_index("seqnum")["ADM0_A3"].replace(NE_RECODE)
    print("Step 3: Kashmir-region cells by NE country:")
    print(kas_admin.value_counts(dropna=False), flush=True)

    # -------------------------------------------------------------------------
    # 6. STEP 4 - ZONAL POPULATION PER HEXAGON (parallel processes; cached per tile)
    # -------------------------------------------------------------------------
    todo = [f for f in shp_tiles if RERUN_ZONAL or not (ZONALDIR / f"{f.stem}.parquet").exists()]
    print(f"Step 4: {len(shp_tiles) - len(todo)} tiles cached, {len(todo)} to run "
          f"({N_WORKERS} workers)", flush=True)
    if todo:
        n_wrap = run_zonal(todo, id_field, rbounds)
        print(f"Step 4: done; antimeridian-split cells in this run = {n_wrap:,}", flush=True)

    zonal = dd.read_parquet(
        ZONALDIR, columns=["seqnum", "cell_area_km2", "wrapped", "pop_sum", "valid_px_cov"]
    )

    # -------------------------------------------------------------------------
    # 7. STEP 5 - COMBINE, ATTACH ATTRIBUTES + POPULATION, APPLY NE OVERRIDE, WRITE
    # -------------------------------------------------------------------------
    BASE = ["seqnum", "longitude", "latitude", "gaul0_code", "match_type", "dist_deg", "dist_km"]
    POP = ["cell_area_km2", "pop_sum", "pop_density", "valid_px_cov", "wrapped"]
    KEEP = ["seqnum", "longitude", "latitude"] + COLS + ["match_type", "dist_deg", "dist_km"] + POP

    # dask Series has notnull(), not notna()
    matched = j1[j1["gaul0_code"].notnull()].assign(
        match_type="intersects", dist_deg=0.0, dist_km=0.0)[BASE]
    near_dd = dd.from_pandas(near[BASE], npartitions=1)
    out = dd.concat([matched, near_dd]).merge(gaul_attr, on="gaul0_code", how="left")

    def apply_overrides(pdf):
        pdf = pdf.copy()
        # (a) Kashmir region: Natural Earth de facto country per cell
        m = pdf["seqnum"].map(kas_admin)
        hit = m.notna()
        pdf.loc[hit, "iso3_admin"] = m[hit]
        pdf.loc[hit, "admin_source"] = "naturalearth"
        # (b) remote islands resolved via Natural Earth
        pdf.loc[pdf["match_type"].eq("ne_remote"), "admin_source"] = "naturalearth_remote"
        # (c) nearest GAUL match too far away and not rescued by NE -> no country at all
        far = pdf["match_type"].eq("nearest") & (pdf["dist_km"] > MAX_NEAREST_KM)
        for c in ["gaul0_code", "gaul0_name", "iso3_code", "iso3_admin", "continent"]:
            pdf.loc[far, c] = pd.NA
        pdf.loc[far, "admin_source"] = "beyond_max_dist"
        return pdf

    out = out.map_partitions(apply_overrides)
    out = out.merge(zonal, on="seqnum", how="left")
    # Outside the raster's latitude extent -> NaN (no data). Inside: nodata = 0 people if
    # NODATA_AS_ZERO, else NaN where the cell has no valid pixels at all.
    inside = (out["latitude"] > rbounds[1]) & (out["latitude"] < rbounds[3])
    keep_mask = inside if NODATA_AS_ZERO else (inside & (out["valid_px_cov"] > 0))
    out["pop_sum"] = out["pop_sum"].where(keep_mask)
    out["pop_density"] = out["pop_sum"] / out["cell_area_km2"]
    out = out[KEEP]

    print(f"Step 5: writing {OUT} ...", flush=True)
    out.to_parquet(OUT, write_index=False, overwrite=True)

    # -------------------------------------------------------------------------
    # 8. QA
    # -------------------------------------------------------------------------
    res = dd.read_parquet(OUT)

    print("\nadmin_source counts:")
    print(res["admin_source"].value_counts().compute())

    print(f"\n(admin_source 'beyond_max_dist' = > {MAX_NEAREST_KM:g} km from GAUL and not found in NE)")
    print("\nmatch_type counts:")
    print(res["match_type"].value_counts().compute())
    rem = res[res["match_type"] == "ne_remote"][["iso3_code", "gaul0_name"]].compute()
    print("\nNE-rescued remote-island cells by country:")
    print(rem.value_counts().head(25))
    far_un = res[res["admin_source"] == "beyond_max_dist"][["seqnum", "latitude", "longitude",
                                                            "dist_km"]].compute()
    far_w = res[res["admin_source"] == "beyond_max_dist"][["wrapped", "pop_sum"]].compute()
    print(f"\nunassigned beyond {MAX_NEAREST_KM:g} km: {len(far_un):,} cells "
          f"({far_w['wrapped'].sum():,} antimeridian-wrapped; total pop_sum "
          f"{far_w['pop_sum'].sum():,.0f}); by rounded location:")
    print(far_un.groupby([far_un.latitude.round(0), far_un.longitude.round(0)])
          .agg(n=("seqnum", "size"), max_km=("dist_km", "max"))
          .sort_values("n", ascending=False).head(20))
    print("\nunassigned iso3_admin by gaul0_code (expect only 100, 110, 270; "
          "beyond_max_dist cells have no gaul0_code):")
    print(res.loc[res["iso3_admin"].isna(), "gaul0_code"].value_counts().compute())

    print("\nnearest-match geodesic distances (km) for cells kept with a GAUL country:")
    print(res.loc[(res["match_type"] == "nearest") & res["iso3_code"].notnull(), "dist_km"]
          .describe().compute())

    print(f"\ncell_area_km2 summary (reference hexagon = {CELL_AREA_KM2:.4f}; pentagons ~5/6):")
    print(res["cell_area_km2"].describe().compute())
    print("centroids with no matching polygon (no zonal result):",
          f"{res['cell_area_km2'].isna().sum().compute():,}")
    print(f"antimeridian-split cells: {res['wrapped'].sum().compute():,}")

    print("\npop_density (people/km2) summary:")
    print(res["pop_density"].describe().compute())
    print(f"pop_density NaN (outside raster extent{'' if NODATA_AS_ZERO else ' or no valid pixels'}): "
          f"{res['pop_density'].isna().sum().compute():,}")
    print("cells with zero valid-pixel coverage (inside extent, set to 0): "
          f"{((res['valid_px_cov'] == 0) & res['pop_density'].notnull()).sum().compute():,}")

    tot_cells = res["pop_sum"].sum().compute()
    tot_rast = raster_total(WORLDPOP)
    print(f"\nsum of cell pop_sum = {tot_cells:,.0f}")
    print(f"WorldPop raster sum = {tot_rast:,.0f}  "
          f"({100 * tot_cells / tot_rast:.2f}% captured by land cells)")

    n_rows = len(res)
    n_unique = res["seqnum"].nunique().compute()
    print(f"\nrows = {n_rows:,}, unique seqnum = {n_unique:,}, duplicates = {n_rows - n_unique:,}")

    # -------------------------------------------------------------------------
    # 9. STEP 6 - GEE JOIN, LANDSCAPE CLASSIFICATION, SUMMARY TABLES
    # -------------------------------------------------------------------------
    run_step6(iso3_to_name)


if __name__ == "__main__":
    main()
