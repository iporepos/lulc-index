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
        "ac": { "name": "Anthropic Cover", "complement_of": "nc" }
      }
    }

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

.. warning::
    The current implementation loads the full source raster into memory.
    For large areas (e.g. full Brazil at 30m), a chunked/windowed approach
    using ``gdal.Warp`` or ``rioxarray`` + ``dask`` is recommended.

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
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform

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

    return spec


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

    for year, src_path in source_files:
        for alias, idx in resolved_indices.items():
            for res_label, res_value in resolutions.items():

                out_path = build_output_path(
                    output, spec_id, alias, year, res_label
                )

                log.info(
                    f"[{completed + skipped + failed + 1}/{total}] "
                    f"{alias} | {year} | {res_label} → {out_path.name}"
                )

                if dry_run:
                    log.info("  [dry-run] skipping computation")
                    skipped += 1
                    continue

                if out_path.exists():
                    log.info("  Output exists, skipping.")
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