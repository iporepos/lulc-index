# lulc-index

Computes Land Use / Land Cover (LULC) index rasters from classified raster
products (e.g. MapBiomas, MODIS MCD12Q1), driven by a JSON spec file.

A LULC index raster expresses, for each pixel at a target resolution, the
fraction of area belonging to a user-defined cluster of land cover classes.
Values range from 0 to 1 and represent continuous spatial gradients of cover
type — natural cover fraction, forest cover fraction, urban fraction, and so on.

---

## How it works

The pipeline performs three operations per output raster:

1. **Binary reclassification** — pixels whose class code belongs to the target
   cluster are set to 1, all others to 0, nodata pixels to NaN.
2. **Spatial aggregation** — the binary mask is resampled to the target
   resolution using average resampling, yielding a 0–1 fraction surface.
3. **Scaling and encoding** — the fraction is multiplied by a scale factor and
   cast to the output dtype for compact storage (e.g. `uint8` with scale factor
   100 stores values 0–100, where 100 = 100% cover).

Two processing modes are available (see [`processing` block](#processing-block-optional)):

- **Conventional** (default) — loads the full source raster into memory, one
  output at a time. Fine for regional rasters (e.g. a state at 30 m).
- **Tiled** — processes the raster in geographic tiles with bounded memory,
  optionally in parallel. Required for very large rasters (e.g. full Brazil at
  30 m, ~22.6 G pixels). Produces pixel-identical results to the conventional
  mode.

---

## Dependencies

```bash
pip install rasterio numpy
```

---

## Spec file

All parameters are declared in a JSON spec file. See [`mapbiomas_col10_example.json`](mapbiomas_col10_example.json) for a working example.

### Top-level structure

| Field | Description |
|---|---|
| `meta` | Provenance metadata: `spec_id`, `description`, `source_legend`, `notes` |
| `source` | Input raster location and file discovery parameters |
| `output` | Output directory, naming template, dtype, scale factor, nodata |
| `resolutions` | Dict of `{ label: degree_value }` target resolutions |
| `indices` | Dict of `{ alias: index_definition }` LULC indices to compute |
| `processing` | *Optional.* Enables tiled processing (`tile_size`, `workers`) |

### `source` block

Source files are resolved in this order of precedence:

| Mode | Fields | Description |
|---|---|---|
| Single file | `dir` + `file` | One static raster, no year dimension |
| Explicit years | `dir` + `pattern` + `years` | Pattern with `{year}` token, filtered to year list |
| Auto-discover | `dir` + `pattern` | Glob all matching files, extract years from filenames |

```json
"source": {
  "dir": "data/mapbiomas",
  "pattern": "mapbiomas_col10_rs_{year}.tif",
  "years": [2018, 2019, 2020, 2021, 2022, 2023, 2024],
  "nodata": 0,
  "crs": "EPSG:4326"
}
```

### `output` block

```json
"output": {
  "dir": "outputs/fractions",
  "naming": "{spec_id}__{alias}__{year}__{resolution}.tif",
  "dtype": "uint8",
  "scale_factor": 100,
  "nodata": 255,
  "crs": "EPSG:4326",
  "qml": "/home/ipo/Data/inputs/default.qml"
}
```

`qml` is optional; see [QGIS style sidecars](#qgis-style-sidecars-optional).

Supported dtypes and their defaults:

| `dtype` | Default `scale_factor` | Valid range | Default `nodata` |
|---|---|---|---|
| `uint8` | 100 | 0–100 | 255 |
| `uint16` | 1000 | 0–1000 | 65535 |
| `float32` | 1 | 0.0–1.0 | -9999.0 |

### `resolutions` block

Keys are labels used in output filenames; values are the target resolution in
the CRS units (degrees for EPSG:4326).

```json
"resolutions": {
  "1ha": 0.001,
  "1km": 0.01
}
```

> **Note on accuracy:** Resolution values are in degrees and are nominal
> approximations referenced to the equator (0.001° ≈ 111 m, 0.01° ≈ 1.11 km).
> At 30°S, pixel area is approximately 7% larger than nominal due to
> latitude-dependent longitude compression. Area-sensitive applications should
> apply a per-pixel latitude correction downstream.

### `indices` block

Indices are keyed by their alias, which appears in output filenames. Each index
either declares a list of class codes (`classes`) or references another index
as its complement (`complement_of`).

```json
"indices": {
  "nc": {
    "name": "Natural Cover",
    "classes": [3, 4, 5, 6, 49, 10, 11, 12, 32, 29, 50, 26, 33, 23]
  },
  "ac": {
    "name": "Anthropic Cover",
    "complement_of": "nc"
  },
  "fc": {
    "name": "Forest Cover",
    "classes": [3, 4, 5, 6, 49, 9]
  }
}
```

A complement index computes `1 - fraction` of its source index. Complement
indices are always resolved after all direct class indices, so the reference is
guaranteed to exist.

### QGIS style sidecars (optional)

A `qml` key may be set at the index level and/or in the `output` block:

```json
"output":  { "...": "...", "qml": "styles/default.qml" },
"indices": {
  "NAT": { "name": "Natural Cover", "classes": [3, 4, 5], "qml": "styles/nat.qml" },
  "URB": { "name": "Urban Area",    "classes": [24] }
}
```

- The index-level `qml` wins; `output.qml` is the fallback default. If both are
  absent, no sidecar is written (silently).
- The QML file is copied next to every output raster with the same basename
  (`BRAZIL_LULC-MB_NAT_1km_1985.tif` → `BRAZIL_LULC-MB_NAT_1km_1985.qml`),
  which QGIS loads automatically when the raster is opened.
- Relative paths are resolved against the spec file directory. A declared file
  that does not exist stops the run when the spec is loaded.
- If an output already exists and is skipped, a missing sidecar is still added.
- One style per index serves all resolutions: values are area fractions,
  independent of pixel size. Class breaks must match the **scaled** values,
  i.e. 0–100 for `uint8` with `scale_factor` 100.

### `processing` block (optional)

When absent, the conventional in-memory mode runs exactly as before (results of
earlier runs stay reproducible). When present, tiled mode is used:

```json
"processing": {
  "tile_size": 100,
  "workers": 4
}
```

| Field | Default | Description |
|---|---|---|
| `tile_size` | 100 | Tile edge in pixels of the **coarsest** resolution in the spec |
| `workers` | 1 | Number of processes computing tiles (only the main process writes) |

**Meaning of `tile_size`.** With resolutions `1ha` = 0.001° and `1km` = 0.01°,
the coarsest is 0.01°, so `tile_size: 100` means 1° × 1° tiles:

| `tile_size` | Ground size | 1km pixels | 1ha pixels | Source pixels (~30 m) | Approx. peak RAM / worker | Tiles (Brazil) |
|---|---|---|---|---|---|---|
| 50 | 0.5° × 0.5° | 50 × 50 | 500 × 500 | ~3.4 M | ~80 MB | ~6,600 |
| 100 | 1° × 1° | 100 × 100 | 1,000 × 1,000 | ~13.8 M | ~300 MB | ~1,650 |
| 200 | 2° × 2° | 200 × 200 | 2,000 × 2,000 | ~55 M | ~1.2 GB | ~420 |

Peak memory scales roughly with `workers` × RAM per worker. Smaller tiles add
overhead; larger tiles use more memory. With source rasters on an external
drive, read throughput may cap the benefit of many workers.

**How it works.**

- All output grids are defined exactly as in the conventional mode (same
  origin, transform and size per resolution).
- Resolutions are fully independent — they need not be multiples of each
  other (e.g. 0.003° and 0.01° work together). At every resolution, an output
  pixel belongs to the tile containing its top-left corner, so tiles never
  overlap or leave seams.
- For each tile, the source window covering all owned output pixels (+1 px
  margin) is read **once** and feeds every index and resolution of that year.
- Tiles holding only source nodata (ocean, neighbouring countries) are written
  as nodata without computation.
- Outputs are written to `*.partial.tif` and renamed only on success, so a
  killed run never leaves a file that a rerun would mistake as finished.

---

## Output files

Output filenames are resolved from the `naming` template:

```
{spec_id}__{alias}__{year}__{resolution}.tif
→ mapbiomas_col10_rs_v1__nc__2023__1km.tif
```

Each output GeoTIFF has provenance embedded as tag metadata:

| Tag | Content |
|---|---|
| `spec_id` | Spec identifier |
| `index_alias` | Index alias |
| `index_name` | Index human-readable name |
| `year` | Source year (absent in single-file mode) |
| `resolution` | Resolution label |
| `scale_factor` | Scale factor applied |
| `created` | ISO 8601 timestamp |
| `processing` | `tiled` (tiled mode only) |
| `tile_size` | Tile size used (tiled mode only) |

If a style is configured, a `.qml` sidecar with the same basename is written
next to each raster.

---

## Usage

```bash
# process all years, indices, and resolutions in the spec
python lulc_index.py -s mapbiomas_col10_example.json

# preview planned outputs without writing any files
python lulc_index.py -s mapbiomas_col10_example.json --dry-run

# process a single index
python lulc_index.py -s mapbiomas_col10_example.json --index NAT

# process a single index and year
python lulc_index.py -s mapbiomas_col10_example.json --index NAT --year 2023

# reprocess existing outputs
python lulc_index.py -s mapbiomas_col10_example.json --overwrite

# full Brazil (spec with a processing block → tiled mode)
python lulc_index.py -s mapbiomas-brazil/mapbiomas_brazil_col10.json
```

### All CLI options

| Flag | Short | Description |
|---|---|---|
| `--spec` | `-s` | Path to the JSON spec file (required) |
| `--index` | `-i` | Process only this index alias |
| `--year` | `-y` | Process only this year |
| `--dry-run` | | Log planned outputs without writing |
| `--overwrite` | | Reprocess existing output files |

---

## Reproducibility

The spec file is the reproducibility contract. Integer class codes in `indices`
are immutable once a spec has been used to produce outputs. Any change to class
membership must result in a new `spec_id` and a new spec file — never an edit
in place. The `spec_id` embedded in every output filename and GeoTIFF tag
ensures full traceability from any raster back to its exact class definitions.

---

## Roadmap

- [x] Windowed/chunked processing for large-area rasters (full Brazil at 30 m) — `processing` block
- [x] QGIS style sidecars — `qml` key
- [ ] MODIS MCD12Q1 example spec
- [ ] Integration into the [`plans`](https://github.com/ipo-exe/plans) hydrology library

---

## License

GNU