#!/usr/bin/env python
"""Download WorldPop Global2 R2025A 100 m CONSTRAINED population rasters (all countries, one
year) and mosaic them into a single global VRT for the DGG zonal-population step.

WorldPop publishes R2025A 100 m data only per country; global mosaics exist only at 1 km.
File pattern (per country, per year):
  https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A/{year}/{ISO3}/v1/100m/
      constrained/{iso3}_pop_{year}_CN_100m_R2025A_v1.tif
If that pattern 404s for a country, the file URL is looked up via the hub REST API (?iso3=).

VRT build (the critical part): country rasters have overlapping bounding boxes, so every
source's nodata MUST be declared transparent (-srcnodata), or later files overwrite earlier ones
with fill values. The script reads the sources' nodata/dtype, passes -srcnodata/-vrtnodata
explicitly, uses THIS environment's gdalbuildvrt with THIS environment's PROJ database, and
validates the result (CRS, nodata, populated control windows) before declaring success.

Paths (year, output folder, VRT name) come from config.py.

Run from the repository root:
    python -u python/04_download_worldpop.py 2>&1 | tee download_worldpop.log
"""

import json
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
import rasterio
from rasterio.windows import from_bounds

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

YEAR = config.POP_YEAR
OUTDIR = Path(config.WORLDPOP_DIR)
VRT = Path(config.WORLDPOP)
GAUL = Path(config.GAUL_GPKG)  # fallback country list
BASE = "https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A"
REST = "https://hub.worldpop.org/rest/data/pop/G2_CN_POP_R25A_100m"
N_THREADS = 6
RETRIES = 3
SKIP_DOWNLOAD = False  # True: only (re)build + validate the VRT from files already on disk
NODATA_OVERRIDE = None  # set (e.g. -99999.0) ONLY if the sources declare no nodata value
CONTROL = {"Gangetic plain": (80, 25, 84, 27), "Nile delta": (30.5, 30.3, 31.5, 31.0)}  # populated
HDR = {"User-Agent": "python"}


# =============================================================================
# DOWNLOAD
# =============================================================================
def url_for(iso3):
    return f"{BASE}/{YEAR}/{iso3.upper()}/v1/100m/constrained/{iso3.lower()}_pop_{YEAR}_CN_100m_R2025A_v1.tif"


def rest_json(url):
    with urlopen(Request(url, headers=HDR), timeout=60) as r:
        return json.load(r)


def list_iso3():
    try:
        iso = sorted({str(d["iso3"]).upper() for d in rest_json(REST).get("data", []) if d.get("iso3")})
        if iso:
            print(f"country list: {len(iso)} ISO3 codes from WorldPop REST API")
            return iso
    except Exception as e:  # noqa: BLE001
        print(f"REST listing failed ({type(e).__name__}: {e}); falling back to GAUL codes")
    import geopandas as gpd
    g = gpd.read_file(GAUL, columns=["iso3_code"], ignore_geometry=True)
    iso = sorted({c.upper() for c in g["iso3_code"].dropna() if not c.startswith("x")})
    print(f"country list: {len(iso)} ISO3 codes from GAUL")
    return iso


def rest_url(iso3):
    """Actual file URL for YEAR via the hub REST API (fallback when the pattern 404s)."""
    try:
        for d in rest_json(f"{REST}?iso3={iso3}").get("data", []):
            if str(d.get("popyear")) == str(YEAR):
                for f in d.get("files", []):
                    if str(f).lower().endswith(".tif"):
                        return f
    except Exception:  # noqa: BLE001
        pass
    return None


def remote_size(url):
    with urlopen(Request(url, method="HEAD", headers=HDR), timeout=60) as r:
        return int(r.headers.get("Content-Length", -1))


def fetch(iso3):
    url, err = url_for(iso3), None
    for attempt in range(1, RETRIES + 1):
        try:
            try:
                size = remote_size(url)
            except HTTPError as e:
                if e.code != 404:
                    raise
                alt = rest_url(iso3)
                if not alt:
                    return iso3, "missing", "not found (404, no REST entry)", 0
                url, size = alt, remote_size(alt)
            dst = OUTDIR / Path(url).name
            if dst.exists() and (size < 0 or dst.stat().st_size == size):
                return iso3, "ok", "cached", dst.stat().st_size
            part = dst.with_name(dst.name + ".part")
            with urlopen(Request(url, headers=HDR), timeout=300) as r, open(part, "wb") as f:
                shutil.copyfileobj(r, f, length=8 << 20)
            if size > 0 and part.stat().st_size != size:
                raise IOError(f"size mismatch {part.stat().st_size} != {size}")
            part.replace(dst)
            return iso3, "ok", f"downloaded ({dst.name})", dst.stat().st_size
        except (HTTPError, URLError, IOError, TimeoutError) as e:
            err = e
            time.sleep(5 * attempt)
    return iso3, "failed", f"FAILED: {type(err).__name__}: {err}", 0


# =============================================================================
# VRT
# =============================================================================
def env_tool(name):
    """Prefer the tool shipped in THIS Python environment (matching PROJ), else PATH."""
    p = Path(sys.prefix) / "bin" / name
    return str(p) if p.exists() else shutil.which(name)


def source_nodata(tifs):
    """The common nodata value declared by all sources (NaN-aware)."""
    vals, dts = [], {}
    for t in tifs:
        with rasterio.open(t) as s:
            vals.append(s.nodata)
            dts[s.dtypes[0]] = dts.get(s.dtypes[0], 0) + 1
    summary = {}
    for v in vals:
        k = "None" if v is None else ("nan" if math.isnan(v) else repr(v))
        summary[k] = summary.get(k, 0) + 1
    print(f"source dtypes: {dts}\nsource nodata: {summary}")
    if len(summary) == 1 and "None" not in summary:
        return vals[0]
    if NODATA_OVERRIDE is not None:
        print(f"using NODATA_OVERRIDE = {NODATA_OVERRIDE}")
        return NODATA_OVERRIDE
    raise SystemExit("sources disagree on / lack a nodata value: inspect one file "
                     "(gdalinfo -stats <file>) and set NODATA_OVERRIDE")


def build_vrt(tifs, nodata):
    exe = env_tool("gdalbuildvrt")
    if not exe:
        raise SystemExit("gdalbuildvrt not found - run: mamba install -c conda-forge gdal")
    env = dict(os.environ)
    proj = Path(sys.prefix) / "share" / "proj"
    if (proj / "proj.db").exists():  # THIS env's PROJ database (avoids version-mismatch errors)
        env["PROJ_DATA"] = env["PROJ_LIB"] = str(proj)
    nd = "nan" if math.isnan(nodata) else repr(float(nodata))
    lst = OUTDIR / "vrt_inputs.txt"
    lst.write_text("\n".join(str(t) for t in tifs))
    cmd = [exe, "-overwrite", "-srcnodata", nd, "-vrtnodata", nd, "-input_file_list", str(lst), str(VRT)]
    print("running:", " ".join(cmd), f"\n  PROJ_DATA={env.get('PROJ_DATA')}", flush=True)
    subprocess.run(cmd, check=True, env=env)


def validate_vrt():
    ok = True
    with rasterio.open(VRT) as s:
        print(f"\nVRT: {s.width}x{s.height}, crs={s.crs}, nodata={s.nodata}, bounds={s.bounds}")
        ok &= s.crs is not None and s.crs.to_epsg() == 4326
        ok &= s.nodata is not None
        for name, b in CONTROL.items():
            a = s.read(1, window=from_bounds(*b, s.transform))
            nd = np.isnan(a) if (s.nodata is None or math.isnan(s.nodata)) else (a == s.nodata)
            v = a[~nd]
            print(f"  {name:15s} >0={np.mean(a > 0):.3f}  nodata={nd.mean():.3f}  "
                  f"negative(valid)={np.mean(v < 0) if v.size else 0:.3f}  people={v[v > 0].sum():,.0f}")
            ok &= bool(np.mean(a > 0) > 0.3) and not bool((v < 0).any())
    print("VRT VALIDATION:", "PASSED" if ok else "FAILED - do not use this VRT")
    return ok


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if not SKIP_DOWNLOAD:
        iso = list_iso3()
        t0, bad = time.time(), {"missing": [], "failed": []}
        with ThreadPoolExecutor(N_THREADS) as ex:
            for i, fut in enumerate(as_completed([ex.submit(fetch, c) for c in iso]), 1):
                c, kind, msg, size = fut.result()
                if kind in bad:
                    bad[kind].append(c)
                print(f"[{i}/{len(iso)}] {c}: {msg} ({size / 1e6:,.1f} MB) | "
                      f"{(time.time() - t0) / 60:.1f} min", flush=True)
        print(f"\nnot on server: {bad['missing']} | FAILED (re-run to retry): {bad['failed']}")
    tifs = sorted(OUTDIR.glob(f"*_{YEAR}_CN_100m_R2025A_v1.tif"))
    print(f"{len(tifs)} country files on disk ({sum(t.stat().st_size for t in tifs) / 1e9:.1f} GB)")
    build_vrt(tifs, source_nodata(tifs))
    if not validate_vrt():
        sys.exit(1)
    print(f"\nwrote {VRT}")


if __name__ == "__main__":
    main()
