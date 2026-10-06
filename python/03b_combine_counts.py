"""
03b_combine_counts.py - combine the Earth Engine count files into one table.

Reads every file matching GEE_RAW_PATTERN in GEE_RAW_DIR (the CSVs exported by
gee/03_extract_counts.js; Parquet files are accepted too), checks that every cell id is unique
(a cell repeated across files must have identical counts), and writes GEE_INPUT: a gzip CSV,
sorted by cell id, with the counts as integers.

Run from the repository root:
    python -u python/03b_combine_counts.py
    python -u python/03b_combine_counts.py path/to/file_or_folder_or_glob ...   (other inputs)
"""

import glob
import gzip
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

ID = "seqnum16"
BANDS = ["n_total", "n_land", "n_freshwater", "n_ice", "n_barren", "n_built", "n_crop", "n_pasture",
         "n_other_intensive_treecrop", "n_nonhabitat", "n_habitat"]
CHUNK_ROWS = 2_000_000


def input_files(args):
    if not args:
        args = [str(Path(config.GEE_RAW_DIR) / config.GEE_RAW_PATTERN)]
    files = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            files += sorted(p.glob(config.GEE_RAW_PATTERN))
        elif any(ch in a for ch in "*?["):
            files += sorted(Path(f) for f in glob.glob(a))
        else:
            files.append(p)
    files = [f for f in files if f.is_file()]
    if not files:
        raise SystemExit(f"no input files found for {args}")
    return files


def read_one(f):
    cols = [ID] + BANDS
    if f.suffix == ".parquet":
        import pyarrow.parquet as pq
        names = pq.read_schema(f).names
        d = pd.read_parquet(f, columns=[c for c in cols if c in names])
    else:
        names = pd.read_csv(f, nrows=0).columns
        d = pd.read_csv(f, usecols=[c for c in cols if c in names], dtype={ID: "string"})
    missing = [c for c in cols if c not in d.columns]
    if missing:
        raise SystemExit(f"{f.name}: missing columns {missing}")
    x = pd.to_numeric(d[ID]).to_numpy(dtype="float64")
    if not (np.isfinite(x).all() and (x == np.floor(x)).all()):
        raise SystemExit(f"{f.name}: non-integer {ID}")
    d[ID] = x.astype("int64")
    for c in BANDS:
        v = pd.to_numeric(d[c]).to_numpy(dtype="float64")
        if not (np.isfinite(v).all() and (v == np.floor(v)).all() and (v >= 0).all()):
            raise SystemExit(f"{f.name}: {c} has missing, negative or non-integer values")
        d[c] = v.astype("int64")
    return d[cols]


def main():
    t0 = time.time()
    files = input_files(sys.argv[1:])
    out = Path(config.GEE_INPUT)
    print(f"{len(files)} input file(s) -> {out}", flush=True)
    d = pd.concat([read_one(f) for f in files], ignore_index=True)
    dup = d[ID].duplicated(keep=False)
    if dup.any():
        conflicting = int((d.loc[dup].groupby(ID)[BANDS].nunique().max(axis=1) > 1).sum())
        if conflicting:
            raise SystemExit(f"{conflicting:,} cells repeated with different counts")
        d = d.drop_duplicates(ID)
    d = d.sort_values(ID, kind="mergesort", ignore_index=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.name}.tmp")
    with gzip.open(tmp, "wt", compresslevel=6, newline="") as fh:
        for i, start in enumerate(range(0, len(d), CHUNK_ROWS)):
            d.iloc[start:start + CHUNK_ROWS].to_csv(fh, index=False, header=(i == 0))
    os.replace(tmp, out)
    print(f"wrote {out}: {len(d):,} cells, {out.stat().st_size / 1e9:.2f} GB "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
