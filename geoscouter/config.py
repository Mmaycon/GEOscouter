"""Paths and cache file names used by the app."""

from __future__ import annotations

import os
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_WORK_DIR = _REPO_ROOT / "data" / "output"


def resolve_work_dir() -> Path:
    """Work directory for uploads, scrape caches, and metadata (override with GEOSCOUTER_WORK_DIR)."""
    override = os.environ.get("GEOSCOUTER_WORK_DIR", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return _DEFAULT_WORK_DIR.resolve()


def ensure_work_dir(path: Path | None = None) -> Path:
    target = path if path is not None else resolve_work_dir()
    try:
        target.mkdir(parents=True, exist_ok=True)
    except PermissionError as exc:
        raise PermissionError(
            f"Cannot create GEOscouter work directory {target}. "
            "Set GEOSCOUTER_WORK_DIR to a writable path (for example ./data/output)."
        ) from exc
    if not os.access(target, os.W_OK):
        raise PermissionError(
            f"GEOscouter work directory is not writable: {target}. "
            "Set GEOSCOUTER_WORK_DIR to a directory you own."
        )
    return target


WORK_DIR = ensure_work_dir()

CACHE_FILES = [
    "geo_webscrap.csv",
    "gds_processed.csv",
    "dubious_gses.csv",
    "filtered_geo_webscrap.csv",
    "curated_gse_catalog.csv",
    "metadata_GSE.xlsx",
    "metadata_filtered_by_word.xlsx",
    "platform_cache.json",
]

GDS_INPUT_NAME = "gds_result.txt"
