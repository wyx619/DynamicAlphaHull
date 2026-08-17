# DynamicAlphaHull

[中文文档](README_zh.md)

`DynamicAlphaHull` constructs geographic range polygons from cleaned WGS 84
occurrence coordinates. It builds a Delaunay-based alpha hull, automatically
selects an alpha that meets coverage and connectivity constraints, and returns
a `terra::SpatVector` polygon ready for spatial analysis.

It is intended for users who need a range polygon that can retain concavities,
gaps, and disjunct occurrences better than a minimum convex hull, without
maintaining legacy scripts or downloading coastline data at runtime.

## Background

This package modernises two older R workflows:

- the Delaunay, alpha-shape, and alpha-hull geometry historically provided by
  `alphahull`;
- the adaptive alpha selection used by `rangeBuilder`, which chooses a range
  from point coverage and the number of polygon parts.

The required functionality was not merely wrapped. It was reorganised into an
independent R package with a modern spatial stack:

- `terra` handles vector geometry, projection, buffering, and overlays;
- no code in this package uses `sf`;
- it does not depend on `alphahull`, `rangeBuilder`, `rnaturalearth`, or a
  network service;
- a pre-dissolved Natural Earth 1:50m land layer is stored as lazy internal
  package data for offline coastline clipping;
- the alpha-hull geometry and adaptive selection semantics are retained, and
  the Delaunay triangulation is sourced from `interp` exactly as in `alphahull`
  so that output matches the reference implementation bit-for-bit.

`DynamicAlphaHull` is deliberately focused: it constructs alpha-hull ranges
from already cleaned occurrence records. It is not a general GIS framework or
an occurrence-record cleaning package.

## Equivalence and performance

The geometry is ported so that output matches the reference `alphahull` /
`rangeBuilder` implementation exactly. In particular, `delvor()` triangulates
with `interp::tri.mesh` — the same single-precision triangulation `alphahull`
uses. This is deliberate: the triangulation edge order determines the alpha-hull
arc order, so using the same triangulation is required for polygon equivalence.
A `deldir`-based triangulation produced the same edge set in a different order,
which changed the resulting polygons; it is therefore no longer used by
`delvor()`.

Performance work is retained elsewhere in the geometry path without changing the
result:

- `ahull()` rejects impossible circle pairs from centre distances before entering
  the original arc-clipping state machine;
- `ashape()` uses vectorised operations and `data.table` aggregation in place
  of per-edge ranking and temporary tables;
- `getDynamicAlphaHull()` re-triangulates on every alpha attempt, mirroring
  `rangeBuilder` exactly so the same `interp` jitter sequence is consumed under
  a fixed seed;
- `getDynamicAlphaHull()` computes the great-circle distance matrix for
  near-duplicate removal once and maintains it incrementally. Removing a point
  does not change the distances among the survivors, so each pass only discards
  the dropped row and column instead of rebuilding the whole O(n^2) matrix. The
  closest-pair selection is unchanged, so the same points are dropped;
- `ah2terra()` ports `rangeBuilder::ah2sf` arc reordering and ring construction,
  and `dummycoor()` uses `interp::in.convex.hull` exactly as `alphahull` does.

End-to-end `getDynamicAlphaHull()` timings on real occurrence records (R 4.6,
Windows, median of warm runs under `set.seed(1)`, `clipToCoast = "terrestrial"`):

| Species | Points | 1.4.0 (s) | 1.5.0 (s) | Speed-up |
| --- | ---: | ---: | ---: | ---: |
| Rubus idaeus | 3,995 | 37.89 | 4.97 | 7.6x |
| Potentilla nivea | 1,654 | 7.73 | 5.00 | 1.5x |
| Neillia incisa | 1,146 | 1.44 | 0.52 | 2.8x |
| Ceanothus fendleri | 801 | 0.75 | 0.21 | 3.6x |

Potentilla nivea needs a 30-attempt alpha search and Rubus idaeus a 15-attempt
search; each attempt re-triangulates in `interp` (the same unoptimised C
triangulation the reference uses) and the sequence is fixed by RNG parity, which
caps how far the largest species can be accelerated. On a fixed alpha the pure-R
geometry path is about 32-50% faster (`ahull()` on the same cleaned records,
e.g. Rubus 0.31 to 0.21 s).

Across a 2,296-species validation against the original rangeBuilder baseline,
2,292 species (99.8%) match with a median area difference of 0.001%. The 1.5.0
optimizations change no output: dropped points, alpha sequences, and polygons
are bit-for-bit identical to 1.4.0. Timings depend on point geometry, the alpha
sequence, and hardware.

## Ellipsoid geometry and the four non-matching species

This package measures distance and area with `terra`'s GeographicLib ellipsoid
and tests point-in-polygon in the same planar space in which the alpha hull is
built. The reference `rangeBuilder` baseline, by contrast, ran with `sf`'s S2
spherical geometry enabled, which models the Earth as a sphere. The ellipsoid
is the more accurate model (the Earth bulges at the equator), so where the two
differ, this package's answer is the more precise one.

Four of 2,296 species diverge from the frozen baseline for library-level rather
than algorithmic reasons:

- `Pouzolzia guatemalana var. nivea` — S2 spherical versus planar
  point-in-polygon in the coverage check;
- `Brosimum rubescens` and `Ficus lutea` — S2/LWGEOM geodesic versus this
  package's great-circle duplicate-point removal, which differ at the
  sub-millimetre level on near-coincident records;
- `Potentilla elegans` — a hull crossing the antimeridian, whose
  self-intersection the `sf` and `terra` builds of GEOS classify differently.

`Ficus lutea` also crashes the reference implementation (an unhandled
`try-error` passed to `sf`), while this package returns a minimum convex hull.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("wyx619/DynamicAlphaHull")
```

The package imports `data.table`, `interp`, `purrr`, and `terra`; missing R
dependencies are installed automatically. On some platforms, installing
`terra` also requires GDAL/PROJ system libraries.

## Quick start

The minimum input is a data frame or matrix with longitude and latitude.
The following creates a range polygon and clips it to land:

```r
library(DynamicAlphaHull)

occurrences <- data.frame(
  Longitude = c(116.0, 116.5, 116.4, 116.1),
  Latitude = c(39.0, 39.0, 39.4, 39.3)
)

range <- getDynamicAlphaHull(
  occurrences,
  buff = 1000,
  clipToCoast = "terrestrial"
)

range$alpha
#> [1] "alpha3"

plot(range$hull, col = "lightblue", border = "steelblue")
points(occurrences$Longitude, occurrences$Latitude, pch = 16)
```

`DynamicAlphaHull` registers a base-R `plot()` method for its output. After
`library(DynamicAlphaHull)`, `plot(range$hull)` works directly; there is no
need to attach `terra` just to draw the range.

The return value is a list with:

- `range$hull`: a WGS 84 `terra::SpatVector` polygon;
- `range$alpha`: the selected alpha label, such as `"alpha0.07"`, or
  `"alphaMCH"` when the minimum-convex-hull fallback is used.

## Adaptive alpha selection

`getDynamicAlphaHull()` starts at `initialAlpha` and adds `alphaIncrement`
until all of the following hold:

1. the resulting polygon is valid;
2. the number of polygon parts is no greater than `partCount`;
3. the polygon intersects at least `fraction` of the input points.

If no candidate meets these conditions by `alphaCap`, the function returns a
buffered minimum convex hull (MCH). Collinear input also uses this fallback.

This example shows a typical multi-step search. Set `verbose = TRUE` to print
each attempted alpha.

```r
set.seed(20260811)
occurrences <- data.frame(
  Longitude = runif(500, 115.8, 116.7),
  Latitude = runif(500, 38.8, 39.6)
)

range <- getDynamicAlphaHull(
  occurrences,
  fraction = 0.98,
  partCount = 1,
  buff = 1000,
  initialAlpha = 0.005,
  alphaIncrement = 0.005,
  alphaCap = 0.08,
  clipToCoast = "no",
  verbose = TRUE
)

range$alpha
#> [1] "alpha0.07"
```

### Parameters and units

| Parameter | Meaning |
| --- | --- |
| `fraction` | Minimum point fraction covered by the range, in `(0, 1]`. |
| `partCount` | Maximum number of disjoint polygon parts. |
| `buff` | Final buffer distance in metres. Buffering is done in Equal Earth (EPSG:8857). |
| `initialAlpha`, `alphaIncrement`, `alphaCap` | Alpha search sequence. Alpha is calculated in the input longitude/latitude plane, so its unit is degrees. |
| `clipToCoast` | `"no"` for no clipping, `"terrestrial"` for land only, and `"aquatic"` for ocean only. |

Alpha is neither a metre-based nor a geodesic distance. The package preserves
the planar alpha-hull definition on longitude/latitude coordinates. Results
for very large extents, antimeridian-crossing records, or polar records should
therefore be interpreted carefully. For regional occurrence data, tune alpha
against the spacing of the points.

### Custom columns and matrices

Use `coordHeaders` when longitude and latitude are not the first two columns.
When the input has exactly two columns, they are used directly:

```r
records <- data.frame(
  species = "example_species",
  lon = c(116.0, 116.5, 116.4, 116.1),
  lat = c(39.0, 39.0, 39.4, 39.3)
)

range <- getDynamicAlphaHull(
  records,
  coordHeaders = c("lon", "lat"),
  clipToCoast = "no"
)

matrixRange <- getDynamicAlphaHull(as.matrix(records[, c("lon", "lat")]))
```

Missing, non-finite, and exact duplicate coordinates are removed. The function
errors when fewer than three unique finite coordinates remain.

## Offline coastline clipping

`clipToCoast = "terrestrial"` and `"aquatic"` use the bundled Natural Earth
1:50m land layer. It is pre-dissolved and stored as lazy internal package data,
then restored and cached only on first use in an R session; no network request,
download, or external data write is made.

```r
land <- loadWorldMap()
plot(land, col = "grey85", border = NA)
plot(range$hull, add = TRUE, border = "firebrick", lwd = 2)
```

## Low-level geometry interface

Most users only need `getDynamicAlphaHull()`. For a fixed alpha, inspection of
intermediate geometry, or reuse of the same Delaunay mesh across several alpha
values, use the lower-level interface:

```r
coordinates <- as.matrix(occurrences[, c("Longitude", "Latitude")])

mesh <- delvor(coordinates)
shape <- ashape(mesh, alpha = 0.05)
hull <- ahull(mesh, alpha = 0.05)
polygon <- ah2terra(hull)

plot(polygon)
```

Passing an existing `delvor` object to `ashape()` or `ahull()` reuses the
triangulation and is preferable when exploring multiple alpha values.

## Project layout

```text
DynamicAlphaHull/
├── R/
│   ├── getDynamicAlphaHull.R  # Adaptive range construction and offline clipping
│   ├── delvor.R               # Delaunay mesh
│   ├── ashape.R               # Alpha-shape edge selection
│   ├── ahull.R                # Alpha-hull arc construction and clipping
│   ├── complement.R           # Complement geometry
│   ├── ah2terra.R             # Arc hull to terra polygon conversion
│   └── geometry.R             # Circle, rotation, arc, and polygon primitives
├── data/
│   └── ne_50m_land.rda        # Lazy internal, dissolved Natural Earth land layer
├── tests/testthat/             # Geometry, Delaunay, and dynamic-range tests
├── man/                        # roxygen2-generated help pages
├── DESCRIPTION                 # Package metadata and dependencies
└── NAMESPACE                   # roxygen2-generated exports
```

Every public function and the registered `plot.SpatVector()` method has
roxygen2 documentation and a small runnable example. Package-level roxygen
imports generate `NAMESPACE`, so implementation code calls its imported
dependencies directly rather than scattering `package::function()` calls. The
test suite covers Delaunay construction, foundation geometry, dynamic ranges,
offline clipping, and degenerate input. During development, regenerate
documentation and run:

```r
roxygen2::roxygenise(".")
testthat::test_local(".")
```

## Dependencies and scope

Runtime dependencies are limited to `data.table`, `interp`, `purrr`, and
`terra`. `interp` supplies the underlying Delaunay triangulation (the same one
used by `alphahull`), while `terra` provides all spatial objects and GIS
operations; this package builds the alpha-shape, alpha-hull, and adaptive range
selection above them.

Complete taxonomic checks, coordinate-precision checks, land/sea consistency
checks, and outlier cleaning before calling this package. A range polygon is a
geometric summary of its input records and chosen parameters; it is not, by
itself, a species' realised distribution or suitable habitat.
