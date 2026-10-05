# ---------------------------------------------------------------------------
# 02_validate_dgg.R - validate the DGG written by 01_build_dgg.R
#
#   0. UNIQUENESS          - no cell id appears in more than one layer.
#   1. AREA RECONCILIATION - informative: total cell area vs Natural Earth land.
#   2. NE INTERIOR         - random points inside Natural Earth land must all
#                            fall in an exported cell.
#   3. NE COASTLINE        - points ON the Natural Earth coastline, likewise.
#   4. OSM INTERIOR        - random points inside (a systematic subset of) OSM
#                            land polygons, likewise.
#   5. OSM COASTLINE       - points on the boundaries of that OSM subset.
#   6. ANTIMERIDIAN        - no polygon part in any seam layer spans > 180 deg
#                            of longitude (Earth Engine would ingest such a
#                            polygon as a globe-spanning band).
#   7. SUPERSET            - optional: every cell of a REFERENCE grid (e.g. an
#                            earlier release) is still present.
#
# Run from the repository root:   Rscript R/02_validate_dgg.R
# ---------------------------------------------------------------------------

library(dggridR)
library(sf)


DATA_DIR     <- Sys.getenv("DATA_DIR", unset = "data")
EXPORT_DIR   <- file.path(DATA_DIR, "dgg")
V2_DIR       <- NA         # optional reference grid folder for check 7 (NA = skip)
PREFIX       <- "dgg_isea3h16_land"
LAND_SHP     <- file.path(DATA_DIR, "inputs", "ne_10m_land", "ne_10m_land.shp")
OSM_GPKG     <- file.path(DATA_DIR, "inputs", "land-polygons-split-4326", "land_polygons.gpkg")
OSM_LAYER    <- "land_polygons"
OSM_EVERY    <- 15L        # use every 15th OSM polygon (by fid) for checks 4-5
N_INTERIOR   <- 200000L
N_COASTLINE  <- 200000L
HEX_AREA_KM2 <- 1.18491
MIN_LAT      <- -90        # must match MIN_LAT in 01_build_dgg.R (-90 = global export)
RUN_CHECKS   <- 0:7        # e.g. 4:7 to rerun only the OSM, antimeridian and superset checks

set.seed(42)   # reproducible sampling


dggs <- dgconstruct(area = 1, metric = TRUE, resround = "nearest")


in_domain <- function(f, min_lat = MIN_LAT) {
  # Layers are owned by centroid latitude: tile lat0 >= min_lat <=> its cells are in the domain.
  b <- basename(f)
  if (grepl("_south_cap", b)) return(min_lat <= -90)
  m <- regmatches(b, regexec("_([ns])([0-9]{2})_[ew][0-9]{3}(_pts\\.csv|\\.shp)$", b))[[1]]
  if (length(m) < 3) return(TRUE)
  lat0 <- as.numeric(m[3]) * ifelse(m[2] == "n", 1, -1)
  lat0 >= min_lat
}


read_exported_seqnums <- function(export_dir = EXPORT_DIR, prefix = PREFIX) {
  files <- list.files(export_dir, pattern = paste0("^", prefix, ".*_pts\\.csv$"),
                      full.names = TRUE)
  n_all <- length(files)
  files <- files[vapply(files, in_domain, logical(1))]
  if (length(files) < n_all) {
    cat(sprintf("Domain (centroid lat >= %g): %d of %d layers\n", MIN_LAT, length(files), n_all))
  }
  if (!length(files)) stop("No ", prefix, "_*_pts.csv files in ", export_dir)
  cat(sprintf("Reading %d tile centroid files from %s ...\n", length(files), export_dir))

  has_dt <- requireNamespace("data.table", quietly = TRUE)
  parts  <- vector("list", length(files))
  for (i in seq_along(files)) {
    parts[[i]] <- if (has_dt) {
      data.table::fread(files[i], select = "seqnum")$seqnum
    } else {
      utils::read.csv(files[i], colClasses = c(seqnum = "numeric"))$seqnum
    }
    if (i %% 100 == 0) cat(sprintf("  %d/%d\n", i, length(files)))
  }
  out <- sort(as.integer(round(as.numeric(unlist(parts, use.names = FALSE)))))
  dups <- sum(diff(out) == 0L)
  cat(sprintf("Exported cells: %s | duplicate ids across layers: %s\n\n",
              format(length(out), big.mark = ","), format(dups, big.mark = ",")))
  attr(out, "dups") <- dups
  out
}


in_set <- function(x, sorted_set) {
  # Membership without hashing 124 million keys.
  idx <- findInterval(x, sorted_set)
  idx > 0L & sorted_set[pmax(idx, 1L)] == x
}


check_area <- function(cells, land) {
  previous_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(TRUE)                       # geodesic area and length
  on.exit(sf::sf_use_s2(previous_s2), add = TRUE)

  land_km2 <- as.numeric(sum(sf::st_area(land))) / 1e6
  coast_km <- as.numeric(sum(sf::st_length(sf::st_cast(land, "MULTILINESTRING")))) / 1e3
  cells_km2 <- length(cells) * HEX_AREA_KM2
  excess <- cells_km2 - land_km2

  cat("---- 1. Area reconciliation (informative; adds OSM land + a membership buffer) ----\n")
  cat(sprintf("Natural Earth land (geodesic) : %15.1f km2\n", land_km2))
  cat(sprintf("Cell area total               : %15.1f km2\n", cells_km2))
  cat(sprintf("Excess over NE land           : %15.1f km2  (%+.2f%%)\n",
              excess, 100 * excess / land_km2))
  cat(sprintf("NE coastline length           : %15.1f km\n", coast_km))
  cat(sprintf("Implied fringe width          : %15.3f km\n", excess / coast_km))
  if (excess < 0) cat("  FAIL: cells cover less than Natural Earth land.\n") else cat("  OK\n")
  cat("\n")
}


check_points <- function(label, points, dggs, cells) {
  coords <- sf::st_coordinates(points)
  # Only points safely inside the domain: a point just north of MIN_LAT can sit in a cell whose
  # centroid is just south of it (cells are ~1 km, ~0.01 deg).
  in_dom <- if (MIN_LAT <= -90) rep(TRUE, nrow(coords)) else coords[, 2] >= MIN_LAT + 0.02
  if (any(!in_dom)) {
    cat(sprintf("(%s sample points south of %g excluded)\n",
                format(sum(!in_dom), big.mark = ","), MIN_LAT + 0.02))
  }
  coords <- coords[in_dom, , drop = FALSE]
  seqnums <- dgGEO_to_SEQNUM(dggs, coords[, 1], coords[, 2])$seqnum
  hit <- in_set(as.integer(seqnums), cells)

  cat(sprintf("---- %s ----\n", label))
  cat(sprintf("Points tested         : %s\n", format(length(hit), big.mark = ",")))
  cat(sprintf("Resolving to a cell   : %s (%.4f%%)\n",
              format(sum(hit), big.mark = ","), 100 * mean(hit)))
  cat(sprintf("Missing               : %s\n", format(sum(!hit), big.mark = ",")))

  if (any(!hit)) {
    misses <- coords[!hit, , drop = FALSE]
    cat("  FAIL: sample locations not covered. First few:\n")
    for (i in seq_len(min(5L, nrow(misses)))) {
      cat(sprintf("    lon %9.4f  lat %8.4f  seqnum %d\n",
                  misses[i, 1], misses[i, 2], seqnums[!hit][i]))
    }
  } else {
    cat("  OK: every sampled location is covered.\n")
  }
  cat("\n")
  invisible(all(hit))
}


check_osm <- function(dggs, cells) {
  q <- sprintf("SELECT * FROM %s WHERE fid %% %d = 0", OSM_LAYER, OSM_EVERY)
  osm <- sf::st_make_valid(sf::st_geometry(sf::st_read(OSM_GPKG, query = q, quiet = TRUE)))
  # Sample in PLANAR lon/lat: with a geographic CRS st_sample() goes through s2, which rejects
  # some rings that GEOS considers valid ("Loop 0 is not valid: Edge ... crosses edge ...").
  # Only the sampled coordinates are used, so dropping the CRS here is harmless.
  sf::st_crs(osm) <- NA
  cat(sprintf("OSM subset for checks 4-5: %s polygons (every %dth)\n",
              format(length(osm), big.mark = ","), OSM_EVERY))
  interior <- sf::st_sample(osm, N_INTERIOR, type = "random", exact = FALSE)
  ok4 <- check_points("4. OSM interior completeness", interior, dggs, cells)
  coast <- sf::st_sample(sf::st_cast(sf::st_boundary(osm), "MULTILINESTRING"), N_COASTLINE,
                         type = "random", exact = FALSE)
  ok5 <- check_points("5. OSM coastline completeness", coast, dggs, cells)
  invisible(ok4 && ok5)
}


check_antimeridian <- function(export_dir = EXPORT_DIR, prefix = PREFIX) {
  # Seam layers: tiles at lon0 = -180 or 170, plus the polar caps.
  shps <- list.files(export_dir, full.names = TRUE,
                     pattern = paste0("^", prefix, ".*(_w180|_e170|_cap)\\.shp$"))
  shps <- shps[vapply(shps, in_domain, logical(1))]
  cat(sprintf("---- 6. Antimeridian (%d seam layers) ----\n", length(shps)))
  n_bad <- 0L
  for (f in shps) {
    g <- sf::st_cast(sf::st_geometry(sf::st_read(f, quiet = TRUE)), "MULTIPOLYGON")
    xy <- sf::st_coordinates(g)                      # columns X, Y, L1 (ring), L2 (part), L3 (feature)
    key <- paste(xy[, "L3"], xy[, "L2"])
    span <- tapply(xy[, "X"], key, function(x) diff(range(x)))
    ymax <- tapply(abs(xy[, "Y"]), key, max)
    bad <- sum(span > 180 & ymax < 89.9)             # pole cells legitimately span all longitudes
    if (bad) cat(sprintf("  FAIL %s: %d polygon parts span > 180 deg\n", basename(f), bad))
    n_bad <- n_bad + bad
  }
  if (n_bad == 0) cat("  OK: no polygon part spans > 180 deg of longitude\n")
  cat("\n")
  invisible(n_bad == 0)
}


check_superset <- function(cells, v2_dir = V2_DIR) {
  cat("---- 7. Superset of reference grid ----\n")
  if (is.na(v2_dir) || !dir.exists(v2_dir)) {
    cat("  skipped (V2_DIR not set / not found)\n\n")
    return(invisible(TRUE))
  }
  old  <- read_exported_seqnums(v2_dir, PREFIX)
  lost <- !in_set(old, cells)
  cat(sprintf("reference cells %s | new cells %s | added %s | reference cells MISSING %s\n",
              format(length(old), big.mark = ","), format(length(cells), big.mark = ","),
              format(length(cells) - (length(old) - sum(lost)), big.mark = ","),
              format(sum(lost), big.mark = ",")))
  if (sum(lost)) cat("  FAIL: the new grid must contain every v2 cell.\n") else cat("  OK\n")
  cat("\n")
  invisible(sum(lost) == 0)
}


main <- function() {
  cells <- read_exported_seqnums()
  if (0 %in% RUN_CHECKS) {
    cat("---- 0. Uniqueness ----\n")
    if (attr(cells, "dups") > 0) cat("  FAIL: duplicate cell ids across layers\n\n") else cat("  OK\n\n")
  }

  previous_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(previous_s2), add = TRUE)

  if (any(1:3 %in% RUN_CHECKS)) {
    cat(sprintf("Reading land layer: %s\n", LAND_SHP))
    land <- sf::st_make_valid(sf::st_union(sf::st_geometry(sf::st_read(LAND_SHP, quiet = TRUE))))
    domain <- sf::st_as_sfc(sf::st_bbox(c(xmin = -180, ymin = MIN_LAT, xmax = 180, ymax = 90),
                                        crs = sf::st_crs(land)))
    land <- sf::st_make_valid(sf::st_intersection(land, domain))   # study domain only
    if (1 %in% RUN_CHECKS) check_area(cells, land)
    if (2 %in% RUN_CHECKS) {
      cat(sprintf("Sampling %s interior points...\n", format(N_INTERIOR, big.mark = ",")))
      interior <- sf::st_sample(land, N_INTERIOR, type = "random", exact = FALSE)
      check_points("2. NE interior completeness", interior, dggs, cells)
    }
    if (3 %in% RUN_CHECKS) {
      cat(sprintf("Sampling %s coastline points...\n", format(N_COASTLINE, big.mark = ",")))
      coast <- sf::st_sample(sf::st_cast(land, "MULTILINESTRING"), N_COASTLINE,
                             type = "random", exact = FALSE)
      check_points("3. NE coastline completeness", coast, dggs, cells)
    }
  }

  if (any(4:5 %in% RUN_CHECKS)) check_osm(dggs, cells)
  if (6 %in% RUN_CHECKS) check_antimeridian()
  if (7 %in% RUN_CHECKS) check_superset(cells)

  invisible(NULL)
}


if (!interactive()) main() else cat("Call main() to run the checks.\n")
