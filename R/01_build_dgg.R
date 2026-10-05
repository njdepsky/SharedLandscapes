# ---------------------------------------------------------------------------
# 01_build_dgg.R - global all-land discrete global grid (DGG)
#
# ISEA3H resolution 16 (DGGRID via dggridR): hexagonal cells of 1.18491 km2,
# 430,467,212 cells on Earth. A cell is "land" (exported) if it touches
#
#     union(Natural Earth 10m land, OSM land polygons), buffered outward by
#     MEMBERSHIP_BUFFER_KM (default 200 m)
#
# The result is a GLOBAL layer (Antarctica included) that can be used for any
# application; study domains (here: centroid latitude >= -60) are applied
# downstream. Reference build (OSM snapshot 2026-09-30): 474 layers,
# 126,509,398 cells.
#
# Output (OUT_DIR), one layer per 10 x 10 degree tile plus two polar caps,
# each cell owned by the tile that contains its centroid:
#   <prefix>_<n|s>LL_<e|w>LLL.zip      zipped shapefile, 'seqnum' as STRING
#                                        (9-digit ids overflow shapefile numbers)
#   <prefix>_<n|s>LL_<e|w>LLL_pts.csv  cell centroids: seqnum, longitude, latitude
#   <prefix>_manifest.csv              layers, cell counts, timings
#
# Method notes
#   - Land is read per tile (OSM from an R-tree-indexed GeoPackage built once),
#     simplified (topology-preserving, SIMPLIFY_KM) and buffered by
#     MEMBERSHIP_BUFFER_KM + SIMPLIFY_KM, so every point of the original land
#     stays at least MEMBERSHIP_BUFFER_KM inside the membership land.
#   - Candidate cells come from dgshptogrid() on a generously buffered sampling
#     polygon; membership is an exact any-overlap test against the land.
#   - Cells crossing the antimeridian are written as east/west parts, each
#     strictly on one side of +/-180; the run STOPS if any part spans > 180 deg
#     (Earth Engine would ingest such a polygon as a globe-spanning band).
#   - Coverage, uniqueness and seam handling are checked by 02_validate_dgg.R.
#
# Run from the repository root:   Rscript R/01_build_dgg.R
# Inputs (see README):  data/inputs/ne_10m_land/ne_10m_land.shp
#                       data/inputs/land-polygons-split-4326/land_polygons.shp
# Override the data folder with the environment variable DATA_DIR.
# ---------------------------------------------------------------------------

library(dggridR)
library(sf)

# 9-digit cell ids must never be written in scientific notation (e.g. "1.07e+08"),
# which turns them into floats downstream.
options(scipen = 999)


# ===========================================================================
# CONFIG
# ===========================================================================

# NOTE: a new directory. Pointing this at a populated OUT_DIR with RESUME
# FALSE stops the run; with RESUME TRUE it continues an interrupted one.
DATA_DIR <- Sys.getenv("DATA_DIR", unset = "data")
OUT_DIR  <- file.path(DATA_DIR, "dgg")
RESUME  <- FALSE

# "bbox"      : cells inside BBOX (fast, easy to check)
# "shapefile" : cells covering the extent of SHP_FILE
# "tiles"     : whole Earth, sharded into TILE_STEP-degree zipped shapefiles
# "land"      : whole Earth masked to LAND_SHP, sharded (see notes below)
# "none"      : define functions only, write nothing
EXPORT_MODE <- "land"

LAYER <- "dgg_isea3h16_nyc"

# Used when EXPORT_MODE == "bbox". ~2x2 degrees is ~31,000 cells at res 16.
BBOX <- list(minlat = 40, minlon = -75, maxlat = 42, maxlon = -73)

# Used when EXPORT_MODE == "shapefile".
SHP_FILE <- NA_character_

# Used when EXPORT_MODE %in% c("tiles", "land"). Degrees per tile.
TILE_STEP <- 10

# Used when EXPORT_MODE == "land".
LAND_SHP <- file.path(DATA_DIR, "inputs", "ne_10m_land", "ne_10m_land.shp")
LAND_PREFIX <- "dgg_isea3h16_land"

# OSM land polygons (WGS84, "split" variant) from
#   https://osmdata.openstreetmap.de/data/land-polygons.html  (ODbL, (c) OpenStreetMap contributors)
# Land used for membership = union(Natural Earth 10m, OSM), buffered by MEMBERSHIP_BUFFER_KM.
OSM_LAND_SHP <- file.path(DATA_DIR, "inputs", "land-polygons-split-4326", "land_polygons.shp")
# Record the snapshot you use (osmdata.openstreetmap.de "Last update" tag); it is printed in the log.
OSM_SNAPSHOT <- "osmdata.openstreetmap.de 'Last update: 2026-09-30T03:38', downloaded 2026-09-30"
OSM_GPKG     <- sub("\\.shp$", ".gpkg", OSM_LAND_SHP)   # spatially indexed copy, built on first use
OSM_LAYER    <- "land_polygons"

# Outward buffer on the land union before the any-overlap membership test, in km. 200 m absorbs
# residual shoreline differences between OSM and the 10 m WorldCover/Dynamic World data (tested
# against 10 m land that Natural Earth misses: 98.5-98.7% captured at 0 m, 98.9% at 200 m,
# 99.1% at 500 m; a 500 m buffer mostly adds empty sea cells).
# Applied in a longitude-stretched frame so it is ~equal in km in both directions (exact at each
# tile's mid-latitude, within about +/-20% across most 10-degree tiles).
MEMBERSHIP_BUFFER_KM <- 0.2

# Southern limit of the EXPORT. Default -90: the layer is a portable, global all-land DGG.
# Tiles are owned by centroid, so a higher value (e.g. -60) skips every tile with lat0 < MIN_LAT
# and the south polar cap. Study domains are applied downstream instead (the GEE extraction
# and post-processing use cells with centroid latitude >= -60).
MIN_LAT <- -90

# Performance: OSM coastlines carry millions of vertices per tile (e.g. Baltic/Finnish
# archipelagos); unioning and buffering them raw is extremely slow. Each land piece is
# simplified (topology-preserving, so small islands are kept) with this tolerance and then
# buffered by MEMBERSHIP_BUFFER_KM + SIMPLIFY_KM, so every point of the ORIGINAL land is still
# at least MEMBERSHIP_BUFFER_KM inside the membership land. Buffers use BUFFER_QUAD_SEGS
# segments per quarter circle (max inward error ~2% of the radius, ~5 m at 250 m).
SIMPLIFY_KM      <- 0.05
BUFFER_QUAD_SEGS <- 8L

# The SAMPLING polygon (land buffered by LAND_BUFFER_DEG, used only to seed candidate cells)
# needs no coastline detail: it is buffered by an extra SAMPLE_SIMPLIFY_DEG and then simplified
# by the same amount, so it still contains every point within LAND_BUFFER_DEG of land while its
# vertex count - which drives dgshptogrid's point-in-polygon cost - drops sharply.
SAMPLE_SIMPLIFY_DEG <- 0.01

# Print per-step timings for every tile (sampling polygon / dgshptogrid / membership test).
VERBOSE_TIMING <- TRUE

# Outward buffer, in degrees, applied to land before cell generation. Cells are
# generated from the buffered polygon so that every cell merely TOUCHING land
# gets sampled, then filtered by a real st_intersects() test against unbuffered
# land. dggetres() reports cell spacing of 1.07509 km (~0.0097 deg at the
# equator); 0.03 deg is ~3x that, a safe margin. Interpreted in DEGREES, which
# requires planar mode - see the PLANAR GEOMETRY note further down.
LAND_BUFFER_DEG <- 0.03

# Spacing, in degrees, of the sample points dggridR uses to find cells.
# Resolution 16 cells are ~1.2 km across (~0.011 deg), so the package default
# of 0.1 silently drops nearly every cell. Keep this well under cell spacing.
CELLSIZE <- 0.004

# Latitude band, measured from each pole, handled by the polar cap builder
# instead of by dgshptogrid sampling. 0.10 deg is ~11 km, comfortably wider
# than the ~0.01 deg band where sampling actually fails.
POLAR_CAP_DEG <- 0.10
POLAR_LAT_STEP <- 0.0005     # cap seeding: latitude spacing, degrees
POLAR_LON_STEP <- 0.25       # cap seeding: longitude spacing, degrees

# Any of "shapefile" (zipped, GEE-ready), "csv" (WKT polygons),
# "centroids" (one point per cell - much smaller).
FORMATS <- c("shapefile", "centroids")

EXPECTED_AREA_KM2 <- 1.18491
EARTH_AREA_KM2    <- 510065621.724


# ===========================================================================
# GRID
# ===========================================================================

dggs <- dgconstruct(
  area     = 1,
  metric   = TRUE,
  resround = "nearest"
)


quiet_planar <- function(expr) {
  # Muffle ONLY the notes that are expected consequences of running in planar
  # mode; every other warning still surfaces. Each entry below was traced to
  # its cause before being added - do not extend this list to silence a
  # warning that has not been explained.
  #
  #   "assumes that they are planar"    GEOS predicates on lon/lat with s2
  #                                     off. Intended: see PLANAR GEOMETRY.
  #   "arc_degrees" /
  #   "assumed to be in decimal degrees"
  #                                     Distances read as degrees, which is
  #                                     what LAND_BUFFER_DEG and CELLSIZE are.
  #   "does not correctly buffer"       The buffer is planar, so 0.03 deg is
  #                                     not a constant distance in km. Fine:
  #                                     it only has to be generous enough to
  #                                     guarantee sample points in every
  #                                     coast-touching cell, and the real
  #                                     membership test is the unbuffered
  #                                     st_intersects(). Over-buffering costs
  #                                     nothing but generation time.
  #   "invalid value range for longlat" shift_lon() translates far-side land
  #                                     by +/-360 for the seam fix, so its
  #                                     coordinates sit outside +/-180 until
  #                                     the sampling polygon is clipped back
  #                                     into range. Expected, and transient.
  #   "does not correctly simplify"     st_simplify() tolerances are in degrees
  #                                     by design (SAMPLE_SIMPLIFY_DEG); the
  #                                     sampling polygon is over-buffered by the
  #                                     same amount, so coverage is preserved.
  expected <- paste(
    "assumes that they are planar",
    "arc_degrees",
    "assumed to be in decimal degrees",
    "does not correctly buffer",
    "invalid value range for longlat",
    "does not correctly simplify",
    sep = "|"
  )
  withCallingHandlers(
    expr,
    # sf wraps long messages over two lines ("assumes that they\nare planar"), so
    # whitespace is normalised before matching.
    warning = function(w) {
      if (grepl(expected, gsub("\\s+", " ", conditionMessage(w)))) invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (grepl(expected, gsub("\\s+", " ", conditionMessage(m)))) invokeRestart("muffleMessage")
    }
  )
}


check_dggs <- function(dggs, expected_area_km2 = EXPECTED_AREA_KM2) {
  info <- dggetres(dggs)
  row  <- info[info$res == dggs$res, , drop = FALSE]

  pick <- function(df, name) if (!is.null(df[[name]])) df[[name]][1] else NA_real_
  area    <- pick(row, "area_km")
  spacing <- pick(row, "spacing_km")

  n_cells <- 10 * 3^dggs$res + 2
  derived <- EARTH_AREA_KM2 / n_cells
  if (is.na(area)) area <- derived

  cat("---- DGG specification ----\n")
  cat(sprintf("Topology / projection : %s / %s\n", dggs$topology, dggs$projection))
  cat(sprintf("Resolution            : %d\n", dggs$res))
  cat(sprintf("Cell area (km2)       : %.5f\n", area))
  cat(sprintf("Cell area, derived    : %.5f\n", derived))
  if (!is.na(spacing)) cat(sprintf("Cell spacing (km)     : %.5f\n", spacing))
  cat(sprintf("Cells on Earth        : %s\n", format(n_cells, big.mark = ",")))

  if (abs(area - expected_area_km2) > 1e-4) {
    stop(sprintf(
      "Cell area %.5f km2 does not match the baseline %.5f km2.",
      area, expected_area_km2
    ))
  }
  cat(sprintf("OK: matches HEX_AREA_KM2 = %.5f\n\n", expected_area_km2))
  invisible(list(res = dggs$res, area_km2 = area, n_cells = n_cells))
}


# ===========================================================================
# GEE-COMPATIBLE WRITERS
# ===========================================================================
# Earth Engine ingests Shapefile (.shp/.dbf/.shx/.prj, zipped) or CSV with a
# WKT/GeoJSON geometry column. Everything must be EPSG:4326. Shapefile field
# names truncate to 10 characters, each sidecar file must stay under 2 GB, the
# zip must hold exactly one shapefile set, and filenames must contain no dots
# beyond the extension. A single table asset caps at 100 million features.

seqnum_column <- function(grid) {
  hit <- which(tolower(names(grid)) == "seqnum")
  if (!length(hit)) stop("No seqnum column found in the grid.")
  grid[[hit[1]]]
}


part_lon_spans <- function(g) {
  # Longitude span of each polygon part (outer ring) of a POLYGON / MULTIPOLYGON sfg.
  if (inherits(g, "POLYGON")) return(diff(range(g[[1]][, 1])))
  vapply(g, function(p) diff(range(p[[1]][, 1])), numeric(1))
}


split_one_antimeridian <- function(g) {
  # The construction Earth Engine ingested correctly for the 1,088 re-extracted seam cells:
  # shift longitudes to 0..360, cut at 180, move the western piece back by -360.
  s <- sf::st_shift_longitude(sf::st_sfc(g, crs = 4326))
  east <- quiet_planar(sf::st_intersection(s, lonlat_box(0, 180, -90, 90)))
  west <- quiet_planar(sf::st_intersection(s, lonlat_box(180, 360, -90, 90)))
  parts <- list()
  if (length(east) && !all(sf::st_is_empty(east))) parts[[length(parts) + 1L]] <- east
  if (length(west) && !all(sf::st_is_empty(west))) parts[[length(parts) + 1L]] <- shift_lon(west, -360)
  u <- quiet_planar(sf::st_union(do.call(c, parts)))
  sf::st_cast(u, "MULTIPOLYGON")[[1]]
}


enforce_antimeridian_split <- function(grid) {
  # Rebuild every seam cell (centre within 2 deg of +/-180, not a pole cell) whose geometry spans
  # > 180 deg of longitude, then STOP if any polygon part still spans > 180 deg.
  geo  <- dgSEQNUM_to_GEO(dggs, as.numeric(seqnum_column(grid)))
  cand <- which(abs(geo$lon_deg) > 178 & abs(geo$lat_deg) < 89.9)
  if (!length(cand)) return(grid)
  gl <- lapply(seq_len(nrow(grid)), function(i) sf::st_geometry(grid)[[i]])
  n_fixed <- 0L
  for (i in cand) {
    if (diff(range(sf::st_coordinates(gl[[i]])[, 1])) > 180) {
      gl[[i]] <- split_one_antimeridian(gl[[i]])
      n_fixed <- n_fixed + 1L
    }
  }
  bad <- cand[vapply(cand, function(i) any(part_lon_spans(gl[[i]]) > 180), logical(1))]
  if (length(bad)) {
    stop(sprintf("%d seam cells still have a polygon part spanning > 180 deg (first seqnum %s)",
                 length(bad), seqnum_column(grid)[bad[1]]))
  }
  sf::st_geometry(grid) <- sf::st_sfc(gl, crs = 4326)
  if (n_fixed) {
    cat(sprintf("  antimeridian: %d seam cells rebuilt as east/west parts; no part spans > 180 deg\n",
                n_fixed))
  }
  grid
}


prepare_for_gee <- function(grid) {
  if (is.na(sf::st_crs(grid))) sf::st_crs(grid) <- 4326
  if (sf::st_crs(grid) != sf::st_crs(4326)) grid <- sf::st_transform(grid, 4326)

  # Cells straddling 180 deg otherwise ingest as ribbons across the globe.
  # Confirmed working: seam tiles come back as 2-part MULTIPOLYGONs with cell
  # area still exactly 1.18491 km2.
  grid <- tryCatch(
    quiet_planar(sf::st_wrap_dateline(
      grid,
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=10"),
      quiet   = TRUE
    )),
    error = function(e) {
      warning("st_wrap_dateline failed; enforcing the split cell by cell: ", conditionMessage(e))
      grid
    }
  )
  # Belt and braces: whatever st_wrap_dateline did, rebuild and CHECK every seam cell.
  grid <- enforce_antimeridian_split(grid)

  # 9-digit ids exceed what a DBF numeric field reliably holds. Character has
  # no such ceiling; GEE reads it as a string property.
  sc <- which(tolower(names(grid)) == "seqnum")
  if (length(sc)) grid[[sc[1]]] <- as.character(grid[[sc[1]]])

  geom_col <- attr(grid, "sf_column")
  is_attr  <- names(grid) != geom_col
  names(grid)[is_attr] <- make.unique(substr(names(grid)[is_attr], 1, 10))
  grid
}


shapefile_sidecars <- function(out_dir, layer) {
  files <- list.files(out_dir, full.names = TRUE)
  keep  <- startsWith(basename(files), paste0(layer, ".")) &
    !endsWith(basename(files), ".zip")
  files[keep]
}


write_gee_shapefile <- function(grid, out_dir, layer) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grid <- prepare_for_gee(grid)

  shp <- file.path(out_dir, paste0(layer, ".shp"))
  for (stale in shapefile_sidecars(out_dir, layer)) unlink(stale)
  sf::st_write(grid, shp, delete_dsn = file.exists(shp), quiet = TRUE)

  sidecars <- shapefile_sidecars(out_dir, layer)
  zipfile  <- file.path(out_dir, paste0(layer, ".zip"))
  if (file.exists(zipfile)) unlink(zipfile)
  utils::zip(zipfile, sidecars, flags = "-j9Xq")   # -j strips directory paths

  oversize <- sidecars[file.size(sidecars) > 2e9]
  if (length(oversize)) {
    warning("Over the 2 GB shapefile limit: ",
            paste(basename(oversize), collapse = ", "))
  }
  cat(sprintf("  shapefile : %s  (%s features)\n",
              basename(zipfile), format(nrow(grid), big.mark = ",")))
  invisible(zipfile)
}


write_gee_csv <- function(grid, out_dir, layer) {
  # Point GEE's "geometry column" advanced option (or the manifest) at `wkt`.
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grid <- prepare_for_gee(grid)

  out <- sf::st_drop_geometry(grid)
  out$wkt <- sf::st_as_text(sf::st_geometry(grid))

  path <- file.path(out_dir, paste0(layer, ".csv"))
  utils::write.csv(out, path, row.names = FALSE, quote = TRUE, na = "")
  cat(sprintf("  csv       : %s  (%s features)\n",
              basename(path), format(nrow(out), big.mark = ",")))
  invisible(path)
}


write_gee_centroids_csv <- function(dggs, seqnums, out_dir, layer) {
  # Exact cell centres. These CSVs, not the shapefiles, are the authoritative
  # record of which cells are in the layer - they never pass through GDAL.
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  geo <- dgSEQNUM_to_GEO(dggs, seqnums)
  out <- data.frame(
    seqnum    = seqnums,
    longitude = geo$lon_deg,
    latitude  = geo$lat_deg
  )
  path <- file.path(out_dir, paste0(layer, "_pts.csv"))
  utils::write.csv(out, path, row.names = FALSE)
  cat(sprintf("  centroids : %s  (%s points)\n",
              basename(path), format(nrow(out), big.mark = ",")))
  invisible(path)
}


write_requested <- function(grid, dggs, out_dir, layer, formats = FORMATS) {
  if ("shapefile" %in% formats) write_gee_shapefile(grid, out_dir, layer)
  if ("csv"       %in% formats) write_gee_csv(grid, out_dir, layer)
  if ("centroids" %in% formats) {
    write_gee_centroids_csv(dggs, seqnum_column(grid), out_dir, layer)
  }
  invisible(NULL)
}


# ===========================================================================
# GRID BUILDERS - single extent
# ===========================================================================

export_bbox <- function(dggs, bbox, out_dir, layer,
                        cellsize = CELLSIZE, formats = FORMATS) {
  cat(sprintf("Building cells for bbox lat [%g, %g], lon [%g, %g]\n",
              bbox$minlat, bbox$maxlat, bbox$minlon, bbox$maxlon))
  grid <- dgrectgrid(
    dggs,
    minlat = bbox$minlat, minlon = bbox$minlon,
    maxlat = bbox$maxlat, maxlon = bbox$maxlon,
    cellsize = cellsize
  )
  write_requested(grid, dggs, out_dir, layer, formats)
  invisible(grid)
}


export_shapefile_extent <- function(dggs, shpfname, out_dir, layer,
                                    cellsize = CELLSIZE, formats = FORMATS) {
  if (is.na(shpfname) || !file.exists(shpfname)) {
    stop("SHP_FILE is not set to an existing shapefile.")
  }
  cat(sprintf("Building cells covering %s\n", shpfname))
  grid <- dgshptogrid(dggs, shpfname, cellsize = cellsize)
  write_requested(grid, dggs, out_dir, layer, formats)
  invisible(grid)
}


tile_layer_name <- function(prefix, lat0, lon0) {
  sprintf("%s_%s%02d_%s%03d",
          prefix,
          if (lat0 < 0) "s" else "n", abs(lat0),
          if (lon0 < 0) "w" else "e", abs(lon0))
}


export_global_tiles <- function(dggs, out_dir, prefix = "dgg_isea3h16",
                                tile_step = TILE_STEP, cellsize = CELLSIZE,
                                formats = "shapefile") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lats  <- seq(-90, 90 - tile_step, by = tile_step)
  lons  <- seq(-180, 180 - tile_step, by = tile_step)
  total <- length(lats) * length(lons)
  done  <- 0L

  for (lat0 in lats) {
    for (lon0 in lons) {
      done  <- done + 1L
      layer <- tile_layer_name(prefix, lat0, lon0)
      if (file.exists(file.path(out_dir, paste0(layer, ".zip")))) next

      cat(sprintf("[%d/%d] %s\n", done, total, layer))
      grid <- tryCatch(
        dgrectgrid(dggs,
                   minlat = lat0, minlon = lon0,
                   maxlat = lat0 + tile_step, maxlon = lon0 + tile_step,
                   cellsize = cellsize),
        error = function(e) {
          warning(layer, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(grid) || nrow(grid) == 0) next
      write_requested(grid, dggs, out_dir, layer, formats)
    }
  }
  invisible(NULL)
}


# ===========================================================================
# GLOBAL LAND-MASKED GRID
# ===========================================================================
# Why this is not "generate the globe, then mask": 430,467,212 hexagons will
# not fit in R memory, and dgshptogrid() over all land at once would need
# roughly 770 million sample points at CELLSIZE. So the Earth is walked one
# tile at a time, ocean-only tiles are skipped outright (about two thirds of
# them), and each land tile is generated and filtered independently.
#
# Membership test: cells are GENERATED from land buffered outward by
# LAND_BUFFER_DEG, which guarantees a sample point inside every cell that so
# much as touches the coast, then KEPT by sf::st_intersects() against the
# unbuffered land. That is a true any-overlap test, not a centroid test - a
# hexagon holding a single square metre of land survives.
#
# Tile ownership: each cell belongs to exactly one tile, decided by its exact
# centre from dgSEQNUM_to_GEO(). Land is tested against a tile box widened by
# the buffer, so a coastal cell straddling a tile edge is still tested against
# complete land rather than a truncated slice. Cells within POLAR_CAP_DEG of
# either pole are excluded here and built by the cap builder instead.
#
# PLANAR GEOMETRY, deliberately. export_global_land() turns s2 off for the
# whole masking run and restores it on exit. Two reasons:
#
#   1. With s2 enabled, st_buffer() on lon/lat data takes its distance in
#      METRES. st_buffer(land, 0.03) would buffer by 3 cm, and the coastal
#      sliver cells this design exists to catch would be dropped silently.
#      Planar mode reads 0.03 as degrees, which is what LAND_BUFFER_DEG means.
#   2. s2 rejects degree-tiling geometry at the poles and the antimeridian.
#
# Over 10-degree tiles the planar/spherical difference in st_intersects() is
# far below one cell width, so nothing is lost by the switch.

prepare_osm <- function(shp = OSM_LAND_SHP, gpkg = OSM_GPKG) {
  # One-off: copy the OSM land polygons into a GeoPackage. Its built-in R-tree spatial index
  # makes the per-tile bounding-box reads below fast; a shapefile would be scanned in full.
  if (!file.exists(gpkg)) {
    if (!file.exists(shp)) stop("OSM_LAND_SHP does not exist: ", shp)
    cat(sprintf("Building spatially indexed copy of OSM land polygons (one-off, minutes):\n  %s\n", gpkg))
    tmp <- paste0(gpkg, ".partial.gpkg")         # an interrupted build never masquerades as complete
    if (file.exists(tmp)) unlink(tmp)
    sf::gdal_utils("vectortranslate", shp, tmp, options = c("-f", "GPKG", "-nln", OSM_LAYER))
    file.rename(tmp, gpkg)
  }
  gpkg
}


load_land <- function(land_shp = LAND_SHP) {
  # Returns list(ne = Natural Earth 10m land as one unioned geometry, osm = GeoPackage path).
  # OSM is far too detailed to union globally; it is read per tile by read_osm_box().
  if (!file.exists(land_shp)) {
    stop("LAND_SHP does not exist: ", land_shp)
  }
  cat(sprintf("Reading land layer: %s\n", land_shp))
  land <- sf::st_read(land_shp, quiet = TRUE)

  if (is.na(sf::st_crs(land))) sf::st_crs(land) <- 4326
  if (sf::st_crs(land) != sf::st_crs(4326)) land <- sf::st_transform(land, 4326)

  land <- sf::st_make_valid(land)
  land <- quiet_planar(sf::st_union(sf::st_geometry(land)))
  osm  <- prepare_osm()
  cat(sprintf("Land layers ready: Natural Earth 10m + OSM land polygons (%s)\n", OSM_SNAPSHOT))
  cat(sprintf("Membership buffer: %.3f km\n", MEMBERSHIP_BUFFER_KM))
  invisible(list(ne = land, osm = osm))
}


read_osm_box <- function(gpkg, box) {
  # OSM land polygons whose bounding box intersects `box` (index-assisted read).
  g <- sf::st_read(gpkg, layer = OSM_LAYER, wkt_filter = sf::st_as_text(box), quiet = TRUE)
  if (!nrow(g)) return(NULL)
  geom <- sf::st_geometry(g)
  if (is.na(sf::st_crs(geom))) sf::st_crs(geom) <- 4326
  sf::st_make_valid(geom)
}


land_in_box <- function(land, box) {
  # Natural Earth land clipped to `box` plus the OSM land polygons whose bounding boxes touch it,
  # returned as separate pieces (NOT unioned: they are buffered piecewise and unioned afterwards,
  # which is far cheaper than unioning raw, vertex-dense coastlines).
  pieces <- list()
  ne <- quiet_planar(sf::st_intersection(land$ne, box))
  if (length(ne) && !all(sf::st_is_empty(ne))) {
    # Clipping can leave line/point slivers where the coastline touches the box edge; keep only
    # polygon parts so those slivers are not buffered into spurious land.
    ne <- suppressWarnings(sf::st_collection_extract(ne, "POLYGON"))
    if (length(ne) && !all(sf::st_is_empty(ne))) pieces[[length(pieces) + 1L]] <- ne
  }
  osm <- read_osm_box(land$osm, box)
  if (!is.null(osm) && length(osm)) pieces[[length(pieces) + 1L]] <- osm
  if (!length(pieces)) return(NULL)
  out <- do.call(c, pieces)
  out[!sf::st_is_empty(out)]
}


km_buffer <- function(geom, km, lat_mid, simplify_km = SIMPLIFY_KM) {
  # Membership land from land pieces: in a frame where longitude is stretched by cos(lat_mid)
  # (so a degree of longitude and of latitude are about equally long), simplify each piece by
  # simplify_km, buffer it by km + simplify_km, union everything, and stretch back.
  if (is.null(geom) || !length(geom)) return(NULL)
  c0 <- max(cos(lat_mid * pi / 180), 0.05)
  g <- sf::st_geometry(geom) * diag(c(c0, 1))
  if (simplify_km > 0) {
    g <- quiet_planar(sf::st_simplify(g, preserveTopology = TRUE, dTolerance = simplify_km / 111.32))
  }
  g <- quiet_planar(sf::st_buffer(g, (km + simplify_km) / 111.32, nQuadSegs = BUFFER_QUAD_SEGS))
  g <- quiet_planar(sf::st_union(g))
  g <- g * diag(c(1 / c0, 1))
  sf::st_crs(g) <- 4326
  sf::st_make_valid(g)
}


tile_box <- function(lat0, lon0, tile_step, pad = 0,
                     clamp_lon = TRUE) {
  # Padding must not push the ring past the poles: a corner at -90.03 wraps
  # around and self-intersects. Longitude clamping is optional because the
  # seam logic below deliberately works in an unclamped, translated frame.
  lo_lat <- max(lat0 - pad, -90)
  hi_lat <- min(lat0 + tile_step + pad, 90)
  lo_lon <- lon0 - pad
  hi_lon <- lon0 + tile_step + pad
  if (clamp_lon) {
    lo_lon <- max(lo_lon, -180)
    hi_lon <- min(hi_lon, 180)
  }

  sf::st_sfc(
    sf::st_polygon(list(matrix(
      c(lo_lon, lo_lat,
        hi_lon, lo_lat,
        hi_lon, hi_lat,
        lo_lon, hi_lat,
        lo_lon, lo_lat),
      ncol = 2, byrow = TRUE
    ))),
    crs = 4326
  )
}


shift_lon <- function(geom, dx) {
  # Translate geometry in longitude, preserving CRS. Used to bring land across
  # the antimeridian into a tile's frame of reference.
  out <- sf::st_geometry(geom) + c(dx, 0)
  sf::st_crs(out) <- 4326
  out
}


lonlat_box <- function(lo_lon, hi_lon, lo_lat, hi_lat) {
  # Explicit four-corner box. tile_box() derives its extent from a single
  # tile_step, which is wrong for the far-side seam strips below: those need a
  # narrow longitude span but the tile's FULL latitude span.
  sf::st_sfc(
    sf::st_polygon(list(matrix(
      c(lo_lon, lo_lat,
        hi_lon, lo_lat,
        hi_lon, hi_lat,
        lo_lon, hi_lat,
        lo_lon, lo_lat),
      ncol = 2, byrow = TRUE
    ))),
    crs = 4326
  )
}


land_for_tile <- function(land, lat0, lon0, tile_step, pad,
                          member_km = MEMBERSHIP_BUFFER_KM) {
  # MEMBERSHIP land for this tile: union(Natural Earth, OSM), INCLUDING anything
  # just across the antimeridian translated into this tile's longitude frame,
  # buffered outward by member_km. Land is read with extra padding so its buffer
  # can reach the tile edge.
  #
  # Without this, a cell centred at 179.99E is tested against land clipped at
  # 180 and cannot see land at 179.99W - a gap along the Chukotka, Aleutian
  # and Fijian seams.
  #
  # Measured effect against Natural Earth 10m: zero cells added or lost across
  # all 22 seam tiles that contain land. NE clips its polygons at +/-180, so
  # seam-crossing cells already overlap land in their own frame. Kept because
  # a land source that does NOT clip at the antimeridian would expose the gap.
  lat_mid <- lat0 + tile_step / 2
  pad <- pad + (member_km / 111.32) / max(cos(lat_mid * pi / 180), 0.05)
  lat_lo <- max(lat0 - pad, -90)
  lat_hi <- min(lat0 + tile_step + pad, 90)

  pieces <- list(land_in_box(land, tile_box(lat0, lon0, tile_step, pad)))

  if (lon0 - pad < -180) {                       # tile reaches past -180
    far <- land_in_box(land, lonlat_box(180 - pad, 180, lat_lo, lat_hi))
    if (!is.null(far)) pieces[[length(pieces) + 1L]] <- shift_lon(far, -360)
  }
  if (lon0 + tile_step + pad > 180) {            # tile reaches past +180
    far <- land_in_box(land, lonlat_box(-180, -180 + pad, lat_lo, lat_hi))
    if (!is.null(far)) pieces[[length(pieces) + 1L]] <- shift_lon(far, 360)
  }

  pieces <- Filter(function(p) !is.null(p) && length(p) > 0 && !all(sf::st_is_empty(p)), pieces)
  if (!length(pieces)) return(NULL)

  # Simplify + buffer each piece, then union (see km_buffer).
  km_buffer(do.call(c, pieces), member_km, lat_mid)
}


cells_owned_by_tile <- function(dggs, grid, lat0, lon0, tile_step,
                                polar_cap_deg = POLAR_CAP_DEG) {
  # Exact cell centres decide ownership, so tiles never double-count a cell.
  geo <- dgSEQNUM_to_GEO(dggs, seqnum_column(grid))
  lon <- geo$lon_deg
  lat <- geo$lat_deg
  lon[lon >= 180]  <- lon[lon >= 180] - 360
  lon[lon < -180]  <- lon[lon < -180] + 360

  hi_lat <- lat0 + tile_step
  hi_lon <- lon0 + tile_step
  in_lat <- if (hi_lat >= 90) lat >= lat0 & lat <= hi_lat else lat >= lat0 & lat < hi_lat
  in_lon <- if (hi_lon >= 180) lon >= lon0 & lon <= hi_lon else lon >= lon0 & lon < hi_lon

  # The caps own everything within polar_cap_deg of a pole.
  outside_cap <- abs(lat) < (90 - polar_cap_deg)

  grid[in_lat & in_lon & outside_cap, , drop = FALSE]
}


cells_for_land_tile <- function(dggs, land, lat0, lon0,
                                tile_step = TILE_STEP, cellsize = CELLSIZE,
                                buffer_deg = LAND_BUFFER_DEG, land_here = NULL) {
  # land_here: membership land already built for this tile (the ocean probe), reused to avoid
  # reading and merging OSM twice.
  if (is.null(land_here)) land_here <- land_for_tile(land, lat0, lon0, tile_step, buffer_deg)
  if (is.null(land_here)) return(NULL)

  # Generate from buffered land so every cell touching the coast is sampled.
  # buffer_deg is in degrees here only because s2 is off; under s2 it would be
  # metres and this buffer would effectively vanish. The sampling polygon is
  # over-buffered then simplified (see SAMPLE_SIMPLIFY_DEG): coverage is kept,
  # vertices are not.
  t0 <- Sys.time()
  sample_poly <- quiet_planar(sf::st_buffer(land_here, buffer_deg + SAMPLE_SIMPLIFY_DEG,
                                            nQuadSegs = BUFFER_QUAD_SEGS))
  sample_poly <- quiet_planar(sf::st_simplify(sample_poly, preserveTopology = TRUE,
                                              dTolerance = SAMPLE_SIMPLIFY_DEG))

  # DGGRID needs valid lon/lat, so the sampling polygon is clipped back into
  # range. The translated far-side land above has already done its job: it
  # put land inside the 0.03 deg strip next to the seam, which is enough for
  # a sample point to land in every seam-crossing cell.
  sample_poly <- quiet_planar(sf::st_intersection(
    sample_poly,
    sf::st_sfc(sf::st_polygon(list(matrix(
      c(-180, -90, 180, -90, 180, 90, -180, 90, -180, -90),
      ncol = 2, byrow = TRUE))), crs = 4326)
  ))
  if (length(sample_poly) == 0 || all(sf::st_is_empty(sample_poly))) return(NULL)

  sample_sf <- sf::st_sf(id = seq_along(sample_poly), geometry = sample_poly)

  t1 <- Sys.time()
  grid <- tryCatch(
    quiet_planar(dgshptogrid(dggs, sample_sf, cellsize = cellsize)),
    error = function(e) {
      warning(sprintf("lat %g lon %g: %s", lat0, lon0, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(grid) || nrow(grid) == 0) return(NULL)

  grid <- cells_owned_by_tile(dggs, grid, lat0, lon0, tile_step)
  if (nrow(grid) == 0) return(NULL)

  # Any-overlap test against the membership land, including the translated
  # far-side pieces so seam cells are judged against complete land. The land is
  # split into its polygon parts so the spatial index can prune per part.
  # st_intersects(x, y) prepares x and indexes y: with x = the large land parts
  # (prepared once) and y = the many small cells, each test is cheap. The
  # reverse order re-tests every cell against unprepared, vertex-dense land.
  t2 <- Sys.time()
  land_parts <- sf::st_cast(land_here, "POLYGON")
  hit_idx <- unique(unlist(quiet_planar(
    sf::st_intersects(land_parts, sf::st_geometry(grid))
  ), use.names = FALSE))
  grid <- grid[sort(hit_idx), , drop = FALSE]
  t3 <- Sys.time()
  if (VERBOSE_TIMING) {
    secs <- function(a, b) as.numeric(difftime(b, a, units = "secs"))
    cat(sprintf("  timing: sampling polygon %.0f s | dgshptogrid %.0f s | membership %.0f s (%s land parts)\n",
                secs(t0, t1), secs(t1, t2), secs(t2, t3), format(length(land_parts), big.mark = ",")))
  }
  if (nrow(grid) == 0) return(NULL)
  grid
}


cells_for_polar_cap <- function(dggs, land, pole = c("south", "north"),
                                cap_deg = POLAR_CAP_DEG,
                                lat_step = POLAR_LAT_STEP,
                                lon_step = POLAR_LON_STEP) {
  # dgshptogrid() samples on a degree grid, which fails within ~0.01 deg of a
  # pole: there a 1.075 km hexagon spans tens of degrees of longitude and
  # almost no latitude, so a latitude row can step straight over it. The
  # previous run lost exactly one cell this way, near the South Pole.
  #
  # Cell ids are seeded from points instead, at a spacing chosen for polar
  # geometry rather than equatorial geometry, then filtered the same way as
  # every other cell.
  pole <- match.arg(pole)
  lats <- if (pole == "south") {
    seq(-90, -90 + cap_deg, by = lat_step)
  } else {
    seq(90 - cap_deg, 90, by = lat_step)
  }
  lons <- seq(-180, 180 - lon_step, by = lon_step)

  cat(sprintf("Seeding %s cap from %s points...\n",
              pole, format(length(lats) * length(lons), big.mark = ",")))

  seeds <- expand.grid(lon = lons, lat = lats, KEEP.OUT.ATTRS = FALSE)
  ids <- unique(dgGEO_to_SEQNUM(dggs, seeds$lon, seeds$lat)$seqnum)
  if (!length(ids)) return(NULL)

  grid <- quiet_planar(dgcellstogrid(dggs, ids))
  if (is.null(grid) || nrow(grid) == 0) return(NULL)

  # Keep only cells whose centre is genuinely inside the cap, so the caps and
  # the regular tiles partition cleanly.
  geo <- dgSEQNUM_to_GEO(dggs, seqnum_column(grid))
  keep <- abs(geo$lat_deg) >= (90 - cap_deg)
  grid <- grid[keep, , drop = FALSE]
  if (nrow(grid) == 0) return(NULL)

  # Membership land for the cap band: union(Natural Earth, OSM), buffered.
  band <- cap_deg + 0.1
  cap_box <- if (pole == "south") lonlat_box(-180, 180, -90, -90 + band) else lonlat_box(-180, 180, 90 - band, 90)
  cap_land <- land_in_box(land, cap_box)
  if (is.null(cap_land)) return(NULL)
  cap_land <- km_buffer(cap_land, MEMBERSHIP_BUFFER_KM, if (pole == "south") -89.9 else 89.9)

  cap_parts <- sf::st_cast(cap_land, "POLYGON")
  hit_idx <- unique(unlist(quiet_planar(
    sf::st_intersects(cap_parts, sf::st_geometry(grid))
  ), use.names = FALSE))
  grid <- grid[sort(hit_idx), , drop = FALSE]
  if (nrow(grid) == 0) return(NULL)

  cat(sprintf("  %s cap: %s land cells\n", pole,
              format(nrow(grid), big.mark = ",")))
  grid
}


export_global_land <- function(dggs, out_dir, land_shp = LAND_SHP,
                               prefix = LAND_PREFIX, tile_step = TILE_STEP,
                               cellsize = CELLSIZE,
                               buffer_deg = LAND_BUFFER_DEG,
                               formats = FORMATS, resume = RESUME) {
  # Expect roughly 124.4 million cells (~29% of the globe). This is a long
  # run - hours, not minutes. Completed tiles are skipped when resume = TRUE.
  previous_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)                 # see the PLANAR GEOMETRY note above
  on.exit(sf::sf_use_s2(previous_s2), add = TRUE)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(out_dir, paste0(prefix, "_manifest.csv"))

  if (file.exists(manifest_path) && !resume) {
    stop("A manifest already exists at\n  ", manifest_path,
         "\nRunning with RESUME = FALSE would skip every completed tile and ",
         "produce nothing.\nEither point OUT_DIR at a new directory, or set ",
         "RESUME <- TRUE to continue that run.")
  }

  land <- load_land(land_shp)

  lats  <- seq(-90, 90 - tile_step, by = tile_step)
  lats  <- lats[lats >= MIN_LAT]                 # no-op at the default MIN_LAT = -90
  lons  <- seq(-180, 180 - tile_step, by = tile_step)
  total <- length(lats) * length(lons)

  manifest <- data.frame(
    layer = character(), lat0 = numeric(), lon0 = numeric(),
    n_cells = integer(), seconds = numeric(),
    stringsAsFactors = FALSE
  )
  if (file.exists(manifest_path)) {
    manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    cat(sprintf("Resuming: %d layers already done\n", nrow(manifest)))
  }

  record <- function(manifest, layer, lat0, lon0, n_cells, seconds) {
    manifest <- rbind(manifest, data.frame(
      layer = layer, lat0 = lat0, lon0 = lon0,
      n_cells = n_cells, seconds = round(seconds, 1),
      stringsAsFactors = FALSE
    ))
    utils::write.csv(manifest, manifest_path, row.names = FALSE)
    manifest
  }

  done    <- 0L
  started <- Sys.time()

  # ---- polar caps first: cheap, and they claim their cells before the
  # ---- regular tiles run, so ownership is unambiguous.
  for (pole in c(if (MIN_LAT <= -90) "south", "north")) {   # south cap only for a global export
    layer <- sprintf("%s_%s_cap", prefix, pole)
    if (layer %in% manifest$layer) next

    t0 <- Sys.time()
    cap <- cells_for_polar_cap(dggs, land, pole)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (is.null(cap)) {
      cat(sprintf("  %s cap: no land cells\n", pole))
      next
    }
    write_requested(cap, dggs, out_dir, layer, formats)
    manifest <- record(manifest, layer,
                       if (pole == "south") -90 else 90, NA, nrow(cap), elapsed)
  }

  # ---- regular tiles
  for (lat0 in lats) {
    for (lon0 in lons) {
      done  <- done + 1L
      layer <- tile_layer_name(prefix, lat0, lon0)

      if (layer %in% manifest$layer) next
      if (file.exists(file.path(out_dir, paste0(layer, ".zip")))) next

      # Cheap ocean rejection before any cell generation. Uses the same
      # seam-aware land so a tile whose only land is across the antimeridian
      # is not skipped.
      probe <- land_for_tile(land, lat0, lon0, tile_step, buffer_deg)
      if (is.null(probe)) next

      cat(sprintf("[%d/%d] %s\n", done, total, layer))
      t0   <- Sys.time()
      grid <- cells_for_land_tile(
        dggs, land, lat0, lon0, tile_step, cellsize, buffer_deg, land_here = probe
      )
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (is.null(grid)) {
        cat("  no land cells\n")
        next
      }
      write_requested(grid, dggs, out_dir, layer, formats)
      manifest <- record(manifest, layer, lat0, lon0, nrow(grid), elapsed)

      cat(sprintf("  running total: %s cells, %.1f min elapsed\n",
                  format(sum(manifest$n_cells), big.mark = ","),
                  as.numeric(difftime(Sys.time(), started, units = "mins"))))
    }
  }

  total_cells <- sum(manifest$n_cells)
  cat(sprintf("\nLayers written     : %d\n", nrow(manifest)))
  cat(sprintf("Land cells total   : %s\n", format(total_cells, big.mark = ",")))
  cat(sprintf("Share of globe     : %.2f%%\n",
              100 * total_cells / (10 * 3^dggs$res + 2)))
  cat(sprintf("Land area implied  : %s km2\n",
              format(round(total_cells * EXPECTED_AREA_KM2), big.mark = ",")))
  cat(sprintf("Manifest           : %s\n", manifest_path))
  if (total_cells > 1e8) {
    cat("\nNOTE: over the 100-million-feature limit for ONE Earth Engine\n")
    cat("table asset. Split on lat0 >= 0 for two assets and combine with\n")
    cat("ee.FeatureCollection.merge(), or use the centroid CSVs.\n")
  }
  cat("\nValidate before ingesting: Rscript R/02_validate_dgg.R\n")
  invisible(manifest)
}


# ===========================================================================
# RUN
# ===========================================================================

main <- function() {
  check_dggs(dggs)

  if (identical(EXPORT_MODE, "none")) {
    cat("EXPORT_MODE is \"none\" - nothing written.\n")
    return(invisible(NULL))
  }

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("Output directory: %s\n", OUT_DIR))

  result <- switch(
    EXPORT_MODE,
    bbox      = export_bbox(dggs, BBOX, OUT_DIR, LAYER),
    shapefile = export_shapefile_extent(dggs, SHP_FILE, OUT_DIR, LAYER),
    tiles     = export_global_tiles(dggs, OUT_DIR),
    land      = export_global_land(dggs, OUT_DIR),
    stop("Unknown EXPORT_MODE: ", EXPORT_MODE)
  )

  cat("\nDone. Files in ", OUT_DIR, "\n", sep = "")
  invisible(result)
}


if (!interactive()) {
  main()
} else {
  cat("Interactive session: call main() to build and write the export.\n")
}


# ===========================================================================
# INGEST
# ===========================================================================
# Small files: Assets tab -> New -> Shape files (or CSV file).
# Large or many files: stage in Cloud Storage, then
#
#   gsutil -m cp data/dgg/*.zip gs://[PLACEHOLDER-bucket]/dgg/
#   earthengine upload table \
#     --asset_id=projects/[PLACEHOLDER-project]/assets/dgg_isea3h16_land_n40_w080 \
#     gs://[PLACEHOLDER-bucket]/dgg/dgg_isea3h16_land_n40_w080.zip
#
# For a WKT CSV add: --primary_geometry_column=wkt
# seqnum ingests as a STRING property; ee.Number.parse() recovers the integer.
# A manifest with several `sources` entries ingests many tiles into one asset,
# still subject to the 100-million-feature limit. For the count extraction upload one
# asset per layer whose tile lat0 >= -60 (into one asset folder); see gee/03_extract_counts.js.
