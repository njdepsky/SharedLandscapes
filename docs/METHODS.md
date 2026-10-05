# Methods — NRI Shared Landscapes on a global DGG

Numbers in brackets refer to the reference run (OSM snapshot 2026-09-30, Dynamic World 2024,
WorldPop R2025A 2024).

## Grid

We use the ISEA3H discrete global grid at resolution 16 (DGGRID via the R package dggridR):
equal-area hexagonal cells of 1.18491 km² on the authalic sphere (12 pentagons). A cell belongs
to the land grid if it touches the union of Natural Earth 10m land and OpenStreetMap land
polygons, buffered outward by 200 m. The buffer absorbs residual shoreline differences between
these vector coastlines and the 10 m land-cover data used below: tested against 10 m land that
Natural Earth alone misses, OSM with a 200 m buffer captures 98.9% of it. The grid is global
[126,509,398 cells in 474 layers]; the study domain is every cell whose centroid lies at or
north of 60°S [114,791,301 cells]. Cells crossing the antimeridian are stored as east/west
parts, and grid completeness, uniqueness and seam handling are validated automatically.

## Land-cover pixel counts

For each cell, Google Earth Engine counts 10 m pixels (EPSG:4326) whose centres fall inside the
hexagon. Valid pixels are those not classed as permanent water in ESA WorldCover 2021 and with
at least one Dynamic World observation in 2024. Land excludes persistent freshwater (Dynamic
World water in ≥ 95% of observations) and persistent snow/ice (≥ 99%), each in clusters of at
least 10 pixels. On land, crop, built-up and bare ground are pixels where that Dynamic World
class occurs in more than 40% of observations; pasture is Global Pasture Watch cultivated
grassland; tree crops are plantations (class 2) in the Spatial Database of Planted Trees (SDPT) Version 2.0. Non-habitat pixels are crop, built-up, pasture or
tree crop; habitat pixels are the remaining land pixels.

## Areas

Land and habitat areas are pixel counts multiplied by the WGS84 ellipsoidal area of one pixel at
the cell's centroid latitude. On 2 million fully inland cells, the counted pixel area reproduces
the measured cell area to within ±0.05% (±0.3% at 70–80°N). Unlike a cell-fraction
apportionment, this does not credit coastal cells with the area of their ocean pixels.

## Population

Population is WorldPop R2025A for 2024 (100 m, constrained), mosaicked from 242 country rasters.
Each cell receives the coverage-weighted sum of the population pixels it overlaps (exactextract,
exact fractional pixel coverage); pixels without data (water and unsettled land in the
constrained product) contribute zero. The cells capture 100.00% of the raster's population
[8.125 billion]. Population density is people per km² of land in the cell.

## Classification

A cell is a **Working Landscape** if it contains land, more than 1% of its land is non-habitat,
and its population density is at least 1 person per km² of land. A Working Landscape is a
**Shared Landscape** if at least 20% of its land is habitat. Land outside Working Landscapes is
**Supporting**; Working Landscapes that are not Shared are **Simplified**. Areas are summed by
country, continent and globe [land 127.61 million km²; Working 26.14 million km² (20.5%); Shared
18.87 million km² (14.8% of land; 72.2% of Working Landscapes)]. Cell counts (n_hex) include
cells with land only.

## Country assignment

Cells are assigned to countries by the position of their centroid in FAO GAUL 2024 (level 0).
Centroids outside every GAUL polygon (mostly coastal cells) take the nearest GAUL country within
50 km (geodesic); beyond that, Natural Earth 10m Admin 0 is used (within 10 km), matched by ISO3
code [2,105 cells on small remote islets remain unassigned, with negligible land]. For areas that
GAUL codes as disputed we apply de facto administration; Abyei, Bir Tawil and the Spratly Islands
have no single administrator and are reported as unassigned [12,692 km² of land in total].
Jammu and Kashmir and adjacent India–China border areas are assigned per cell using Natural
Earth's de facto boundaries (India, Pakistan, China; Siachen Glacier to India). Elsewhere GAUL's
country coding is used as published. Continents follow GAUL.
