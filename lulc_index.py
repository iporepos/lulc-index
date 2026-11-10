"""
lulc_index.py
-------------
Computes Land Use / Land Cover (LULC) index rasters from a classified raster
product (e.g. MapBiomas, MODIS MCD12Q1) driven by a JSON spec file.

A LULC index raster expresses, for each pixel at a target resolution, the
fraction of area belonging to a user-defined cluster of land cover classes.
Values range from 0 to 1 (stored as scaled integers for storage efficiency)
and represent continuous spatial gradients of cover type — for example,
natural cover fraction, forest cover fraction, or urban fraction.

The pipeline performs three operations per output:

1. **Binary reclassification** — pixels belonging to the target class cluster
   are set to 1, all others to 0, nodata pixels to NaN.
2. **Spatial aggregation** — the binary mask is resampled to the target
   resolution using average resampling, yielding a 0–1 fraction surface.
3. **Scaling and encoding** — the fraction is multiplied by a scale factor
   and cast to the output dtype for compact storage.

**Spec file structure:**

.. code-block:: json

    {
      "meta":        { "spec_id": "...", "description": "...", "notes": "..." },
      "source":      { "dir": "...", "pattern": "..._{year}.tif", "years": [...] },
      "output":      { "dir": "...", "naming": "...", "dtype": "uint8", ... },
      "resolutions": { "1ha": 0.001, "1km": 0.01 },
      "indices": {
        "nc": { "name": "Natural Cover", "classes": [3, 4, 5] },
        "ac": { "name": "Anthropic Cover", "complement_of": "nc",
                "qml": "styles/ac.qml" }
      }
    }

**QGIS style sidecars (optional):**

Each index may declare a ``qml`` key pointing to a QGIS style file, and the
``output`` block may declare a ``qml`` key used as default for indices that
do not declare their own. When set, the QML is copied next to every output
raster with the same basename (``<output>.qml``), which QGIS loads
automatically when the raster is opened. Relative paths are resolved against
the spec file directory. Style class breaks must match the *scaled* values
(e.g. 0–100 for ``uint8`` with ``scale_factor`` 100).

**Source resolution modes** (in order of precedence):

1. ``file`` key — single static file, no year dimension.
2. ``pattern`` + ``years`` list — explicit year filter.
3. ``pattern`` only — auto-discover all matching files, extract years from filenames.

Output filenames are resolved from the ``naming`` template, e.g.::

    {spec_id}__{alias}__{year}__{resolution}.tif
    → mapbiomas_col10_rs_v1__nc__2023__1km.tif

Provenance (spec_id, alias, year, resolution, scale_factor, timestamp) is
embedded in each output GeoTIFF as tag metadata.

.. note::
    Resolution values in the spec are in degrees (EPSG:4326). Labels such as
    ``1ha`` and ``1km`` are nominal approximations referenced to the equator.
    At 30°S, pixel area is ~7% larger than nominal. Area-sensitive applications
    should apply a per-pixel latitude correction downstream.

**Processing modes:**

- *Conventional* (no ``processing`` block): each output is computed by
  loading the full source raster into memory. Kept as-is for backward
  compatibility and reproducibility of earlier results. Not viable for very
  large rasters (full Brazil at 30 m needs well over 100 GB per output).
- *Tiled* (``processing`` block present)::

      "processing": { "tile_size": 100, "workers": 1 }

  ``tile_size`` is the tile edge in pixels of the **coarsest** resolution of
  the spec (100 × 0.01° = 1° × 1° tiles). Resolutions stay fully independent:
  at every resolution an output pixel belongs to the tile containing its
  top-left corner, so tiles never overlap or leave seams. Each source window
  is read once per tile and feeds all indices and resolutions of that year.
  ``workers`` (default 1) sets the number of processes computing tiles; only
  the main process writes. Output grids are identical to the conventional
  mode. Outputs are written to ``*.partial.tif`` and renamed on success.

.. note::
    This script is intended to be merged into the ``plans`` hydrology library
    once the spec format and pipeline architecture have stabilized.

**Usage:**

.. code-block:: bash

    python lulc_index.py -s specs/mapbiomas_col10_example.json
    python lulc_index.py -s specs/mapbiomas_col10_example.json --dry-run
    python lulc_index.py -s specs/mapbiomas_col10_example.json --index nc
    python lulc_index.py -s specs/mapbiomas_col10_example.json --index nc --year 2023
    python lulc_index.py -s specs/mapbiomas_col10_example.json --overwrite

:dependencies: rasterio, numpy
"""

import argparse
import json
import logging
import math
import multiprocessing as mp
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform
from rasterio.windows import Window, from_bounds
from rasterio.windows import transform as window_transform

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VALID_DTYPES = {"uint8", "uint16", "float32"}
"""Supported output data types."""

DTYPE_DEFAULTS = {
    "uint8":   {"scale_factor": 100,    "nodata": 255},
    "uint16":  {"scale_factor": 1000,   "nodata": 65535},
    "float32": {"scale_factor": 1,      "nodata": -9999.0},
}
"""Default scale factor and nodata sentinel per dtype.

The nodata value is chosen as the maximum representable value for integer
types (outside the valid scaled range) and a conventional sentinel for float.
"""

NUMPY_DTYPES = {
    "uint8":   np.uint8,
    "uint16":  np.uint16,
    "float32": np.float32,
}
"""Mapping from dtype string to numpy dtype."""


# ---------------------------------------------------------------------------
# Spec loading and validation
# ---------------------------------------------------------------------------

def validate_source(source: dict) -> None:
    """Validate the ``source`` block of a spec for contradictions and missing fields.

    :param source: The ``source`` dict from the loaded spec.
    :raises ValueError: If the source block contains contradictory or invalid fields.
    """
    has_file = "file" in source
    has_pattern = "pattern" in source

    if not has_file and not has_pattern:
        raise ValueError("Source must have either 'file' or 'pattern'.")

    if has_file and has_pattern:
        raise ValueError("Source cannot have both 'file' and 'pattern'.")

    if has_file and "years" in source:
        raise ValueError(
            "Source has 'file' (single file mode) — 'years' is not applicable. "
            "Use 'pattern' for multi-year sources."
        )

    if has_pattern:
        pattern = source["pattern"]
        if "{year}" not in pattern:
            raise ValueError(
                f"Pattern '{pattern}' contains no '{{year}}' token. "
                "Use 'file' for a single static source, "
                "or add '{year}' to the pattern for multi-year sources."
            )


def load_spec(spec_path: Path) -> dict:
    """Load and validate a JSON spec file.

    Fills in default ``scale_factor`` and ``nodata`` values for the output
    block if they are not explicitly declared, based on the declared ``dtype``.

    :param spec_path: Path to the JSON spec file.
    :returns: Validated spec dict with output defaults applied.
    :raises ValueError: If required keys are missing or the dtype is unsupported.
    :raises FileNotFoundError: If ``spec_path`` does not exist.
    """
    with open(spec_path) as f:
        spec = json.load(f)

    required_keys = {"source", "output", "resolutions", "indices"}
    missing = required_keys - spec.keys()
    if missing:
        raise ValueError(f"Spec missing required keys: {missing}")

    validate_source(spec["source"])

    dtype = spec["output"].get("dtype", "uint8")
    if dtype not in VALID_DTYPES:
        raise ValueError(f"Invalid dtype '{dtype}'. Must be one of {VALID_DTYPES}.")

    defaults = DTYPE_DEFAULTS[dtype]
    spec["output"].setdefault("scale_factor", defaults["scale_factor"])
    spec["output"].setdefault("nodata", defaults["nodata"])

    _resolve_qml_paths(spec, spec_path.parent)
    _validate_processing(spec)

    return spec


def _validate_processing(spec: dict) -> None:
    """Validate the optional ``processing`` block and apply its defaults.

    Absent block = conventional (in-memory) mode, nothing to check.

    :param spec: Loaded spec dict, modified in place.
    :raises ValueError: If ``tile_size`` or ``workers`` are not positive integers.
    """
    proc = spec.get("processing")
    if proc is None:
        return
    if not isinstance(proc, dict):
        raise ValueError("'processing' must be an object.")
    unknown = set(proc) - {"tile_size", "workers"}
    if unknown:
        log.warning(f"Unknown keys in 'processing' ignored: {sorted(unknown)}")
    proc.setdefault("tile_size", 100)
    proc.setdefault("workers", 1)
    for key in ("tile_size", "workers"):
        val = proc[key]
        if not isinstance(val, int) or isinstance(val, bool) or val < 1:
            raise ValueError(f"processing.{key} must be a positive integer, got {val!r}.")


def _resolve_qml_paths(spec: dict, base_dir: Path) -> None:
    """Resolve and validate optional ``qml`` style paths in place.

    Looks at ``output.qml`` (default for all indices) and ``indices.<alias>.qml``.
    Relative paths are resolved against ``base_dir`` (the spec file directory).

    :param spec: Loaded spec dict, modified in place.
    :param base_dir: Directory used to resolve relative paths.
    :raises FileNotFoundError: If a declared QML file does not exist.
    """
    targets = [("output", spec["output"])]
    targets += [(f"indices.{a}", idx) for a, idx in spec["indices"].items()]
    for where, block in targets:
        if not block.get("qml"):
            continue
        path = Path(block["qml"]).expanduser()
        if not path.is_absolute():
            path = (base_dir / path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"QML style not found ({where}.qml): {path}")
        block["qml"] = path


# ---------------------------------------------------------------------------
# Source file discovery
# ---------------------------------------------------------------------------

def _extract_year_from_filename(path: Path, pattern: str) -> int | None:
    """Extract a year integer from a filename using a pattern with a ``{year}`` token.

    :param path: File path whose name will be parsed.
    :param pattern: Filename pattern containing a ``{year}`` token,
        e.g. ``mapbiomas_col10_rs_{year}.tif``.
    :returns: Extracted year as an integer, or ``None`` if parsing fails.
    """
    prefix, suffix = pattern.split("{year}")
    name = path.name
    if name.startswith(prefix) and name.endswith(suffix):
        year_str = name[len(prefix): len(name) - len(suffix)]
        try:
            return int(year_str)
        except ValueError:
            return None
    return None


def resolve_source_files(source: dict) -> list[tuple[int | None, Path]]:
    """Resolve source raster files from the spec ``source`` block.

    Source resolution hierarchy (first match wins):

    1. ``file`` key — single static file, ``year=None``.
    2. ``pattern`` + ``years`` list — resolve each year explicitly.
    3. ``pattern`` only — glob directory, auto-extract years from filenames.

    :param source: The ``source`` dict from the loaded spec.
    :returns: List of ``(year, path)`` tuples. ``year`` is ``None`` in
        single-file mode.
    :raises FileNotFoundError: If no matching files are found.
    :raises ValueError: If year extraction fails for all matched files.
    """
    src_dir = Path(source["dir"])

    # option 1 — single static file
    if "file" in source:
        path = src_dir / source["file"]
        if not path.exists():
            raise FileNotFoundError(f"Source file not found: {path}")
        log.info(f"Single-file mode: {path}")
        return [(None, path)]

    pattern = source["pattern"]
    years = source.get("years")

    # option 2 — explicit year list
    if years:
        files = []
        for year in years:
            path = src_dir / pattern.format(year=year)
            if not path.exists():
                log.warning(f"Source file not found, skipping: {path}")
                continue
            files.append((year, path))
        if not files:
            raise FileNotFoundError(
                f"No source files found for any of the specified years: {years}"
            )
        return files

    # option 3 — auto-discover
    glob_pattern = pattern.replace("{year}", "*")
    found = sorted(src_dir.glob(glob_pattern))
    if not found:
        raise FileNotFoundError(
            f"No files matched pattern '{glob_pattern}' in '{src_dir}'"
        )

    files = []
    for path in found:
        year = _extract_year_from_filename(path, pattern)
        if year is None:
            log.warning(f"Could not parse year from filename, skipping: {path.name}")
            continue
        files.append((year, path))

    if not files:
        raise ValueError(
            f"Matched files but could not extract year from any filename "
            f"using pattern '{pattern}'."
        )

    log.info(
        f"Auto-discovered {len(files)} source file(s): "
        f"years {[y for y, _ in files]}"
    )
    return files


# ---------------------------------------------------------------------------
# Index resolution
# ---------------------------------------------------------------------------

def resolve_indices(indices: dict) -> dict:
    """Resolve index definitions, expanding complement references.

    Processes direct class indices first (pass 1), then resolves
    ``complement_of`` references (pass 2), ensuring dependencies are always
    available before they are needed.

    :param indices: The ``indices`` dict from the loaded spec, keyed by alias.
    :returns: Resolved indices dict. Each entry contains:

        - ``name`` (str) — human-readable index name.
        - ``classes`` (list[int]) — target class codes.
        - ``is_complement`` (bool) — whether to apply ``1 - fraction``.
        - ``source_alias`` (str, complement only) — alias of the referenced index.
        - ``qml`` (Path | None) — QGIS style file for the sidecar, if declared.

    :raises ValueError: If a ``complement_of`` reference points to an unknown alias
        or to another complement index.
    """
    resolved = {}

    # pass 1 — direct class indices
    for alias, idx in indices.items():
        if "classes" in idx:
            resolved[alias] = {
                "name": idx["name"],
                "classes": idx["classes"],
                "is_complement": False,
                "qml": idx.get("qml"),
            }

    # pass 2 — complement indices
    for alias, idx in indices.items():
        if "complement_of" in idx:
            source_alias = idx["complement_of"]
            if source_alias not in resolved:
                raise ValueError(
                    f"Index '{alias}' references complement_of '{source_alias}' "
                    f"which is not a direct class index."
                )
            resolved[alias] = {
                "name": idx["name"],
                "classes": resolved[source_alias]["classes"],
                "is_complement": True,
                "source_alias": source_alias,
                "qml": idx.get("qml"),
            }

    return resolved


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_binary_mask(data: np.ndarray, classes: list[int]) -> np.ndarray:
    """Compute a binary float32 mask from a classified raster array.

    :param data: 2D integer array of land cover class codes.
    :param classes: List of class codes to include (set to 1).
    :returns: Float32 array of the same shape: ``1.0`` where the pixel class
        is in ``classes``, ``0.0`` otherwise.
    """
    return np.isin(data, classes).astype(np.float32)


def compute_fraction(
    src_path: Path,
    classes: list[int],
    is_complement: bool,
    target_res: float,
    src_nodata: int,
    output: dict,
) -> tuple[np.ndarray, dict]:
    """Compute a LULC fraction raster from a classified source raster.

    Reads the source raster, builds a binary mask for the target classes,
    optionally complements it, propagates nodata as NaN, then resamples to
    the target resolution using average resampling. The result is scaled and
    cast to the output dtype.

    :param src_path: Path to the source classified raster (GeoTIFF).
    :param classes: List of integer class codes defining the index.
    :param is_complement: If ``True``, the fraction is computed as
        ``1 - fraction`` of the target classes (i.e. all other classes).
    :param target_res: Target output resolution in the source CRS units
        (degrees for EPSG:4326).
    :param src_nodata: Nodata value in the source raster.
    :param output: Output configuration dict with keys ``dtype``,
        ``scale_factor``, and ``nodata``.
    :returns: Tuple of ``(result_array, rasterio_profile)``.

    .. warning::
        This function reads the entire source raster into memory. For very
        large rasters, consider a windowed or chunked approach.
    """
    dtype = output["dtype"]
    scale_factor = output["scale_factor"]
    out_nodata = output["nodata"]
    np_dtype = NUMPY_DTYPES[dtype]

    with rasterio.open(src_path) as src:
        data = src.read(1)
        src_nodata_mask = data == src_nodata

        mask = compute_binary_mask(data, classes)

        if is_complement:
            mask = 1.0 - mask

        # nodata pixels must not contribute to the resampled average
        mask[src_nodata_mask] = np.nan

        dst_transform, dst_width, dst_height = calculate_default_transform(
            src.crs, src.crs,
            src.width, src.height,
            *src.bounds,
            resolution=target_res,
        )

        dst = np.full((dst_height, dst_width), np.nan, dtype=np.float32)

        reproject(
            mask, dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=src.crs,
            resampling=Resampling.average,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )

        valid = ~np.isnan(dst)
        result = np.full((dst_height, dst_width), out_nodata, dtype=np_dtype)

        if dtype == "float32":
            result[valid] = dst[valid].astype(np_dtype)
        else:
            result[valid] = np.round(dst[valid] * scale_factor).astype(np_dtype)

        profile = {
            "driver": "GTiff",
            "dtype": dtype,
            "width": dst_width,
            "height": dst_height,
            "count": 1,
            "crs": src.crs,
            "transform": dst_transform,
            "nodata": out_nodata,
            "compress": "lzw",
            "tiled": True,
            "blockxsize": 256,
            "blockysize": 256,
        }

        return result, profile


# ---------------------------------------------------------------------------
# Output path builder
# ---------------------------------------------------------------------------

def build_output_path(
    output: dict,
    spec_id: str,
    alias: str,
    year: int | None,
    resolution_label: str,
) -> Path:
    """Resolve the output file path from the naming template.

    If ``year`` is ``None`` (single-file mode), the ``{year}`` token and any
    adjacent separators (``_`` or ``__``) are stripped so the filename stays
    clean.

    :param output: Output configuration dict containing ``dir`` and ``naming``.
    :param spec_id: Spec identifier string (from ``meta.spec_id``).
    :param alias: Index alias key (e.g. ``'nc'``).
    :param year: Year integer, or ``None`` for single-file mode.
    :param resolution_label: Resolution label key (e.g. ``'1km'``).
    :returns: Resolved output :class:`~pathlib.Path`.
    """
    naming = output["naming"]

    if year is None:
        naming = re.sub(r'_*\{year\}_*', '_', naming).strip('_')
        name = naming.format(spec_id=spec_id, alias=alias, resolution=resolution_label)
    else:
        name = naming.format(
            spec_id=spec_id,
            alias=alias,
            year=year,
            resolution=resolution_label,
        )

    return Path(output["dir"]) / name


# ---------------------------------------------------------------------------
# Tiled processing
# ---------------------------------------------------------------------------

_EPS = 1e-9

_W_SRC = None
"""Source dataset opened once per worker process (tiled mode)."""
_W_CTX = None
"""Shared, read-only job context per worker process (tiled mode)."""


def _encode(dst: np.ndarray, output: dict) -> np.ndarray:
    """Scale a 0–1 float fraction array (NaN = nodata) to the output dtype.

    Same encoding as :func:`compute_fraction`.
    """
    np_dtype = NUMPY_DTYPES[output["dtype"]]
    valid = ~np.isnan(dst)
    result = np.full(dst.shape, output["nodata"], dtype=np_dtype)
    if output["dtype"] == "float32":
        result[valid] = dst[valid].astype(np_dtype)
    else:
        result[valid] = np.round(dst[valid] * output["scale_factor"]).astype(np_dtype)
    return result


def _plan_grids(src, resolutions: dict) -> dict:
    """Full output grid per resolution, identical to the conventional mode.

    :returns: ``{label: (transform, width, height, res)}``.
    """
    grids = {}
    for label, res in resolutions.items():
        t, w, h = calculate_default_transform(
            src.crs, src.crs, src.width, src.height, *src.bounds, resolution=res,
        )
        grids[label] = (t, w, h, res)
    return grids


def _tile_edges(n_tiles: int, tile_ext: float, res: float, size: int) -> list[int]:
    """Pixel offsets where each tile starts along one axis, plus the end.

    A pixel belongs to the tile containing its leading edge, so edges are
    ``ceil(i * tile_ext / res)``; consecutive tiles share edges (no gaps,
    no overlaps). The last edge is clamped to the grid size.
    """
    edges = [min(size, max(0, math.ceil(i * tile_ext / res - _EPS))) for i in range(n_tiles)]
    edges.append(size)
    return edges


def _plan_tiles(src, grids: dict, tile_size: int) -> list[dict]:
    """Split all output grids into common geographic tiles.

    :returns: List of jobs, each with ``src_window`` (or ``None`` when the
        tile is outside the source) and ``out`` = ``{label: Window}``.
    """
    coarsest = max(g[3] for g in grids.values())
    tile_ext = tile_size * coarsest
    n_cols = max(math.ceil(g[1] * g[3] / tile_ext - _EPS) for g in grids.values())
    n_rows = max(math.ceil(g[2] * g[3] / tile_ext - _EPS) for g in grids.values())

    col_edges = {k: _tile_edges(n_cols, tile_ext, g[3], g[1]) for k, g in grids.items()}
    row_edges = {k: _tile_edges(n_rows, tile_ext, g[3], g[2]) for k, g in grids.items()}
    full = Window(0, 0, src.width, src.height)

    jobs = []
    for ti in range(n_rows):
        for tj in range(n_cols):
            out = {}
            xs, ys = [], []
            for k, (t, w, h, res) in grids.items():
                c0, c1 = col_edges[k][tj], col_edges[k][tj + 1]
                r0, r1 = row_edges[k][ti], row_edges[k][ti + 1]
                if c1 <= c0 or r1 <= r0:
                    continue
                win = Window(c0, r0, c1 - c0, r1 - r0)
                out[k] = win
                left, top = t * (c0, r0)
                right, bottom = t * (c1, r1)
                xs += [left, right]
                ys += [top, bottom]
            if not out:
                continue

            # source window covering every owned output pixel, +1 px margin
            sw = from_bounds(min(xs), min(ys), max(xs), max(ys), src.transform)
            c0 = math.floor(sw.col_off) - 1
            r0 = math.floor(sw.row_off) - 1
            c1 = math.ceil(sw.col_off + sw.width) + 1
            r1 = math.ceil(sw.row_off + sw.height) + 1
            sw = Window(c0, r0, c1 - c0, r1 - r0)
            try:
                sw = sw.intersection(full)
                if sw.width <= 0 or sw.height <= 0:
                    sw = None
            except Exception:  # WindowError: no overlap
                sw = None
            jobs.append({"src_window": sw, "out": out})
    return jobs


def _process_tile(job: dict, src=None, ctx=None) -> dict | None:
    """Compute all index × resolution blocks of one tile.

    :returns: ``{(alias, label): uint/float array}`` or ``None`` when the
        tile holds only source nodata (caller writes nodata blocks).
    """
    src = src if src is not None else _W_SRC
    ctx = ctx if ctx is not None else _W_CTX
    sw = job["src_window"]
    if sw is None:
        return None

    data = src.read(1, window=sw)
    nodata_mask = data == ctx["src_nodata"]
    if nodata_mask.all():
        return None
    sw_transform = window_transform(sw, src.transform)

    results = {}
    base_cache = {}
    for alias, idx in ctx["indices"].items():
        key = tuple(sorted(idx["classes"]))
        if key not in base_cache:
            base_cache[key] = compute_binary_mask(data, idx["classes"])
        mask = base_cache[key]
        mask = (1.0 - mask) if idx["is_complement"] else mask.copy()
        mask[nodata_mask] = np.nan

        for label, win in job["out"].items():
            t = ctx["grids"][label][0]
            dst = np.full((int(win.height), int(win.width)), np.nan, dtype=np.float32)
            reproject(
                mask, dst,
                src_transform=sw_transform,
                src_crs=src.crs,
                dst_transform=window_transform(win, t),
                dst_crs=src.crs,
                resampling=Resampling.average,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
            results[(alias, label)] = _encode(dst, ctx["output"])
    return results


def _worker_init(src_path: str, ctx: dict) -> None:
    """Pool initializer: open the source once per worker process."""
    global _W_SRC, _W_CTX
    _W_SRC = rasterio.open(src_path)
    _W_CTX = ctx


def _provenance_tags(spec_id, alias, idx, res_label, output, year, extra=None) -> dict:
    """GeoTIFF tag metadata written to each output."""
    tags = dict(
        spec_id=spec_id,
        index_alias=alias,
        index_name=idx["name"],
        resolution=res_label,
        scale_factor=str(output["scale_factor"]),
        created=datetime.now().isoformat(timespec="seconds"),
    )
    if year is not None:
        tags["year"] = str(year)
    if extra:
        tags.update(extra)
    return tags


def run_tiled_year(
    spec_id: str,
    year: int | None,
    src_path: Path,
    indices: dict,
    resolutions: dict,
    source: dict,
    output: dict,
    processing: dict,
    dry_run: bool,
) -> tuple[int, int, int]:
    """Compute every pending index × resolution output of one source file tile by tile.

    The source is read once per tile for all outputs. Outputs that already
    exist are skipped (style sidecar backfilled), as in the conventional mode.

    :returns: ``(completed, skipped, failed)`` counts for this year.
    """
    tile_size = processing["tile_size"]
    workers = processing["workers"]

    pending = {}
    skipped = 0
    for alias, idx in indices.items():
        for res_label in resolutions:
            out_path = build_output_path(output, spec_id, alias, year, res_label)
            qml_src = idx.get("qml") or output.get("qml")
            if dry_run:
                log.info(f"  [dry-run] {alias} | {year} | {res_label} → {out_path.name}"
                         + (f" (style {qml_src.name})" if qml_src else ""))
                skipped += 1
                continue
            if out_path.exists():
                log.info(f"  {out_path.name} exists, skipping.")
                if qml_src and not out_path.with_suffix(".qml").is_file():
                    write_qml_sidecar(qml_src, out_path)
                    log.info("  Style sidecar added.")
                skipped += 1
                continue
            pending[(alias, res_label)] = (out_path, qml_src)

    with rasterio.open(src_path) as src:
        grids = _plan_grids(src, resolutions)
        jobs = _plan_tiles(src, grids, tile_size)
        n_src = sum(1 for j in jobs if j["src_window"] is not None)
        log.info(
            f"  Tiled: {year} | tile_size {tile_size} "
            f"({tile_size * max(resolutions.values()):g}° tiles) | {len(jobs)} tiles "
            f"({n_src} overlap the source) | workers {workers} | "
            f"{len(pending)} output(s) pending"
        )
        if dry_run or not pending:
            return 0, skipped, 0

        # only compute what is pending
        ctx = {
            "indices": {a: indices[a] for a in {a for a, _ in pending}},
            "grids": grids,
            "src_nodata": source.get("nodata", 0),
            "output": output,
        }
        labels_needed = {lab for _, lab in pending}
        for j in jobs:
            j["out"] = {k: w for k, w in j["out"].items() if k in labels_needed}

        handles = {}
        partials = {}
        try:
            for (alias, label), (out_path, _) in pending.items():
                t, w, h, _res = grids[label]
                profile = {
                    "driver": "GTiff",
                    "dtype": output["dtype"],
                    "width": w,
                    "height": h,
                    "count": 1,
                    "crs": src.crs,
                    "transform": t,
                    "nodata": output["nodata"],
                    "compress": "lzw",
                    "tiled": True,
                    "blockxsize": 256,
                    "blockysize": 256,
                    "BIGTIFF": "IF_SAFER",
                }
                part = out_path.with_name(out_path.stem + ".partial.tif")
                part.unlink(missing_ok=True)
                partials[(alias, label)] = part
                handles[(alias, label)] = rasterio.open(part, "w", **profile)

            np_dtype = NUMPY_DTYPES[output["dtype"]]

            def write(job, res):
                for (alias, label), dst in handles.items():
                    win = job["out"].get(label)
                    if win is None:
                        continue
                    if res is None:
                        arr = np.full((int(win.height), int(win.width)),
                                      output["nodata"], dtype=np_dtype)
                    else:
                        arr = res[(alias, label)]
                    dst.write(arr, 1, window=win)

            step = max(1, len(jobs) // 20)
            if workers == 1:
                stream = ((j, _process_tile(j, src, ctx)) for j in jobs)
                pool = None
            else:
                pool = mp.get_context("fork").Pool(
                    workers, initializer=_worker_init, initargs=(str(src_path), ctx),
                )
                # ordered results keep writes in row-major order (GTiff block cache friendly)
                stream = zip(jobs, pool.imap(_process_tile, jobs, chunksize=1))
            try:
                for n, (job, res) in enumerate(stream, 1):
                    write(job, res)
                    if n % step == 0 or n == len(jobs):
                        log.info(f"  tiles {n}/{len(jobs)} ({100 * n / len(jobs):.0f}%)")
            finally:
                if pool is not None:
                    pool.close()
                    pool.join()

            for (alias, label), dst in handles.items():
                dst.update_tags(**_provenance_tags(
                    spec_id, alias, indices[alias], label, output, year,
                    extra={"processing": "tiled", "tile_size": str(tile_size)},
                ))
                dst.close()
        except Exception as e:
            log.error(f"  FAILED ({year}): {e}")
            for dst in handles.values():
                if not dst.closed:
                    dst.close()
            for part in partials.values():
                part.unlink(missing_ok=True)
            return 0, skipped, len(pending)

    completed = 0
    for key, (out_path, qml_src) in pending.items():
        partials[key].replace(out_path)
        log.info(f"  Written: {out_path}")
        if write_qml_sidecar(qml_src, out_path):
            log.info(f"  Style:   {out_path.with_suffix('.qml').name}")
        completed += 1
    return completed, skipped, 0


def write_qml_sidecar(qml_src: Path | None, raster_path: Path) -> Path | None:
    """Copy a QGIS style file next to a raster as ``<raster stem>.qml``.

    :param qml_src: Template QML path, or ``None`` to do nothing.
    :param raster_path: Output raster path the sidecar belongs to.
    :returns: Path of the written sidecar, or ``None`` if no style was set.
    """
    if qml_src is None:
        return None
    dst = raster_path.with_suffix(".qml")
    shutil.copyfile(qml_src, dst)
    return dst


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def run_pipeline(
    spec: dict,
    filter_index: str | None = None,
    filter_year: int | None = None,
    dry_run: bool = False,
) -> None:
    """Execute the LULC index pipeline for all year × index × resolution combinations.

    Iterates over all combinations defined by the spec (or filtered subsets),
    computes each fraction raster, and writes it to the output directory with
    provenance embedded as GeoTIFF tag metadata.

    Existing output files are skipped unless ``--overwrite`` is passed at the
    CLI level (handled in :func:`main`).

    :param spec: Validated spec dict as returned by :func:`load_spec`.
    :param filter_index: If given, process only this index alias.
    :param filter_year: If given, process only this year.
    :param dry_run: If ``True``, log planned outputs without writing any files.
    """
    spec_id = spec.get("meta", {}).get("spec_id", "unknown_spec")
    source = spec["source"]
    output = spec["output"]
    resolutions = spec["resolutions"]

    source_files = resolve_source_files(source)
    resolved_indices = resolve_indices(spec["indices"])

    if filter_year:
        source_files = [(y, p) for y, p in source_files if y == filter_year]
    if filter_index:
        if filter_index not in resolved_indices:
            raise ValueError(
                f"Index '{filter_index}' not found in spec. "
                f"Available: {list(resolved_indices.keys())}"
            )
        resolved_indices = {filter_index: resolved_indices[filter_index]}

    if not source_files:
        log.warning("No source files to process after filtering.")
        return

    out_dir = Path(output["dir"])
    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    total = len(source_files) * len(resolved_indices) * len(resolutions)
    log.info(
        f"Pipeline: {len(source_files)} year(s) × "
        f"{len(resolved_indices)} index/indices × "
        f"{len(resolutions)} resolution(s) = {total} output(s)"
    )

    completed = 0
    skipped = 0
    failed = 0
    processing = spec.get("processing")

    for year, src_path in source_files:
        if processing is not None:
            c, s_, f = run_tiled_year(
                spec_id, year, src_path, resolved_indices, resolutions,
                source, output, processing, dry_run,
            )
            completed += c
            skipped += s_
            failed += f
            continue

        for alias, idx in resolved_indices.items():
            for res_label, res_value in resolutions.items():

                out_path = build_output_path(
                    output, spec_id, alias, year, res_label
                )
                qml_src = idx.get("qml") or output.get("qml")

                log.info(
                    f"[{completed + skipped + failed + 1}/{total}] "
                    f"{alias} | {year} | {res_label} → {out_path.name}"
                )

                if dry_run:
                    if qml_src:
                        log.info(f"  [dry-run] style sidecar from {qml_src}")
                    log.info("  [dry-run] skipping computation")
                    skipped += 1
                    continue

                if out_path.exists():
                    log.info("  Output exists, skipping.")
                    # backfill the style sidecar for outputs from earlier runs
                    if qml_src and not out_path.with_suffix(".qml").is_file():
                        write_qml_sidecar(qml_src, out_path)
                        log.info("  Style sidecar added.")
                    skipped += 1
                    continue

                try:
                    result, profile = compute_fraction(
                        src_path=src_path,
                        classes=idx["classes"],
                        is_complement=idx["is_complement"],
                        target_res=res_value,
                        src_nodata=source.get("nodata", 0),
                        output=output,
                    )

                    with rasterio.open(out_path, "w", **profile) as dst:
                        dst.write(result, 1)
                        tags = dict(
                            spec_id=spec_id,
                            index_alias=alias,
                            index_name=idx["name"],
                            resolution=res_label,
                            scale_factor=str(output["scale_factor"]),
                            created=datetime.now().isoformat(timespec="seconds"),
                        )
                        if year is not None:
                            tags["year"] = str(year)
                        dst.update_tags(**tags)

                    log.info(f"  Written: {out_path}")
                    if write_qml_sidecar(qml_src, out_path):
                        log.info(f"  Style:   {out_path.with_suffix('.qml').name}")
                    completed += 1

                except Exception as e:
                    log.error(f"  FAILED: {e}")
                    failed += 1

    log.info(
        f"Done. Completed: {completed} | Skipped: {skipped} | Failed: {failed}"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    :returns: Parsed :class:`~argparse.Namespace` with attributes
        ``spec``, ``index``, ``year``, ``dry_run``, and ``overwrite``.
    """
    parser = argparse.ArgumentParser(
        description="Compute LULC index rasters from a JSON spec file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python lulc_index.py -s specs/mapbiomas_col10_example.json
  python lulc_index.py -s specs/mapbiomas_col10_example.json --dry-run
  python lulc_index.py -s specs/mapbiomas_col10_example.json --index nc
  python lulc_index.py -s specs/mapbiomas_col10_example.json --index nc --year 2023
  python lulc_index.py -s specs/mapbiomas_col10_example.json --overwrite
        """,
    )
    parser.add_argument(
        "-s", "--spec",
        required=True,
        type=Path,
        help="Path to the JSON spec file.",
    )
    parser.add_argument(
        "-i", "--index",
        default=None,
        help="Process only this index alias (e.g. 'nc'). Default: all indices.",
    )
    parser.add_argument(
        "-y", "--year",
        default=None,
        type=int,
        help="Process only this year. Default: all years in spec.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned outputs without writing any files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files. Default: skip existing.",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for the CLI.

    Loads and validates the spec file, then delegates to :func:`run_pipeline`.
    When ``--overwrite`` is set, patches :meth:`~pathlib.Path.exists` to
    bypass the skip-existing logic inside the pipeline loop.
    """
    args = parse_args()

    if not args.spec.exists():
        log.error(f"Spec file not found: {args.spec}")
        sys.exit(1)

    log.info(f"Loading spec: {args.spec}")
    spec = load_spec(args.spec)

    spec_id = spec.get("meta", {}).get("spec_id", "unknown")
    log.info(f"Spec ID: {spec_id}")

    if args.overwrite:
        original_exists = Path.exists
        Path.exists = lambda self: False  # type: ignore[method-assign]
        try:
            run_pipeline(spec, args.index, args.year, args.dry_run)
        finally:
            Path.exists = original_exists  # type: ignore[method-assign]
    else:
        run_pipeline(spec, args.index, args.year, args.dry_run)


if __name__ == "__main__":
    main()