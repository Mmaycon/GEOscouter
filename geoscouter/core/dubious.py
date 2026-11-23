"""Flag GSEs whose supplementary files suggest a different platform than the declared GPL."""

from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

from geoscouter.core.curation import unified_platform_name

GENERIC = "Generic/unspecified"
SPATIAL_UNSPECIFIED = "Spatial (unspecified)"

# Ordered: the first signature reaching min_hits wins, so specific techs precede plain 10x.
TECH_FILE_SIGNATURES: list[tuple[str, list[str], int]] = [
    (
        "Xenium",
        [
            r"(?<!detected_)transcripts\.(parquet|csv)",
            r"(^|[_\-.])cells\.(parquet|csv)",
            r"cell_feature_matrix",
            r"nucleus_boundaries",
            r"morphology(_focus|_mip)?(\.ome)?\.tif",
            r"experiment\.xenium",
        ],
        2,
    ),
    (
        "CosMx",
        [
            r"exprmat_file",
            r"metadata_file",
            r"fov_positions",
            r"polygons",
            r"tx_file",
        ],
        2,
    ),
    (
        "MERFISH/MERSCOPE",
        [
            r"cell_by_gene",
            r"detected_transcripts",
            r"cell_metadata",
        ],
        2,
    ),
    (
        "Visium",
        [
            r"scalefactors_json",
            r"tissue_positions",
            r"aligned_fiducials",
            r"tissue_(hires|lowres)_image",
            r"detected_tissue_image",
        ],
        2,
    ),
    (
        "scRNA-seq",
        [
            r"barcodes\.tsv",
            r"(features|genes)\.tsv",
            r"matrix\.mtx",
        ],
        3,
    ),
    ("DNA methylation", [r"\.idat(\.gz)?$"], 1),
    ("ChIP-seq", [r"narrowpeak|broadpeak", r"\.bw(\.gz)?$|bigwig"], 1),
]

_COMPILED = [
    (tech, [re.compile(tok, re.IGNORECASE) for tok in tokens], min_hits)
    for tech, tokens, min_hits in TECH_FILE_SIGNATURES
]

# Techs that have dedicated GPLs (or are expected to be recoverable by platform), so a
# generic sequencer GPL is suspicious. Plain 10x scRNA-seq and ChIP-seq always sit on
# sequencer GPLs and are only flagged when the GPL names a different technology.
GPL_IDENTIFIABLE_TECHS = {"Xenium", "CosMx", "MERFISH/MERSCOPE", "Visium", "DNA methylation"}

SPATIAL_TECHS = {"Xenium", "CosMx", "MERFISH/MERSCOPE", "Visium"}

_DECLARED_MAP = {
    "Xenium": "Xenium",
    "CosMx": "CosMx",
    "Visium": "Visium",
    "MERFISH": "MERFISH/MERSCOPE",
    "MERSCOPE": "MERFISH/MERSCOPE",
    "seqFISH": "seqFISH",
    "Slide-seq": "Slide-seq",
    "scRNA-seq": "scRNA-seq",
    "snRNA-seq": "scRNA-seq",
    "Spatial transcriptomics": SPATIAL_UNSPECIFIED,
}

_METHYLATION_RE = re.compile(r"methylation|infinium|450k|\bepic\b", re.IGNORECASE)

# (declared, file signature) pairs that are consistent even though the labels differ:
# a Visium or unspecified spatial series may only ship its 10x matrix files.
_COMPATIBLE = {
    ("Visium", "scRNA-seq"),
    (SPATIAL_UNSPECIFIED, "scRNA-seq"),
} | {(SPATIAL_UNSPECIFIED, tech) for tech in SPATIAL_TECHS}


def detect_file_technology(filenames: Iterable[str]) -> tuple[str | None, list[str]]:
    """Infer a technology from a GSE's supplementary file names."""
    names = [str(n).strip() for n in filenames if isinstance(n, str) and str(n).strip()]
    for tech, patterns, min_hits in _COMPILED:
        matched = []
        for pat in patterns:
            hit = next((n for n in names if pat.search(n)), None)
            if hit is not None:
                matched.append(pat.pattern)
        if len(matched) >= min_hits:
            return tech, matched
    return None, []


def declared_technology(
    platform_title: str,
    platform_technology: str = "",
    platforms: str = "",
) -> str:
    """Map the declared GPL to the same vocabulary as ``detect_file_technology``."""
    title = str(platform_title or "")
    tech = str(platform_technology or "")
    label = unified_platform_name(title, tech, str(platforms or ""))
    if label in _DECLARED_MAP:
        return _DECLARED_MAP[label]
    if _METHYLATION_RE.search(f"{title} {tech}"):
        return "DNA methylation"
    return GENERIC


def _dubious_reason(declared: str, file_tech: str | None) -> str:
    if file_tech is None or declared == file_tech or (declared, file_tech) in _COMPATIBLE:
        return ""
    if declared == GENERIC:
        return "generic platform" if file_tech in GPL_IDENTIFIABLE_TECHS else ""
    return "conflicting platform"


def _first_text(values: pd.Series) -> str:
    vals = values.dropna().astype(str).str.strip()
    vals = vals[vals != ""]
    return vals.iloc[0] if not vals.empty else ""


def flag_dubious_gses(df: pd.DataFrame) -> pd.DataFrame:
    """One row per Series comparing the declared platform with the file-signature technology."""
    columns = [
        "Series", "Platforms", "Platform_title", "Platform_labels", "Assay_hint",
        "Declared_tech", "File_signature_tech", "Matched_tokens", "N_files",
        "Dubious", "Reason",
    ]
    if df is None or df.empty or "Series" not in df.columns:
        return pd.DataFrame(columns=columns)

    rows = []
    for gse, grp in df.groupby("Series", sort=False):
        text = {
            col: _first_text(grp[col]) if col in grp.columns else ""
            for col in ("Platforms", "Platform_title", "Platform_technology",
                        "Platform_labels", "Assay_hint")
        }
        files = grp.get("Supplementary file", pd.Series(dtype=str)).dropna().astype(str)
        file_tech, matched = detect_file_technology(files.unique())
        declared = declared_technology(
            text["Platform_title"], text["Platform_technology"], text["Platforms"]
        )
        reason = _dubious_reason(declared, file_tech)
        rows.append({
            "Series": gse,
            "Platforms": text["Platforms"],
            "Platform_title": text["Platform_title"],
            "Platform_labels": text["Platform_labels"],
            "Assay_hint": text["Assay_hint"],
            "Declared_tech": declared,
            "File_signature_tech": file_tech or "",
            "Matched_tokens": "; ".join(matched),
            "N_files": int(files.nunique()),
            "Dubious": bool(reason),
            "Reason": reason,
        })
    out = pd.DataFrame(rows, columns=columns)
    return out.sort_values(["Dubious", "Series"], ascending=[False, True]).reset_index(drop=True)
