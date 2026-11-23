"""Tests for platform vs file-signature mismatch detection."""

from __future__ import annotations

import unittest

import pandas as pd

from geoscouter.core.dubious import (
    GENERIC,
    declared_technology,
    detect_file_technology,
    flag_dubious_gses,
)

VISIUM_FILES = [
    "GSM8522424_01_SP1_aligned_fiducials.jpg.gz",
    "GSM8522424_01_SP1_barcodes.tsv.gz",
    "GSM8522424_01_SP1_detected_tissue_image.jpg.gz",
    "GSM8522424_01_SP1_features.tsv.gz",
    "GSM8522424_01_SP1_matrix.mtx.gz",
    "GSM8522424_01_SP1_scalefactors_json.json.gz",
    "GSM8522424_01_SP1_tissue_hires_image.png.gz",
]
XENIUM_FILES = [
    "GSM1_cells.parquet.gz",
    "GSM1_transcripts.parquet.gz",
    "GSM1_cell_feature_matrix.h5",
    "GSM1_morphology.ome.tif.gz",
]
SCRNA_FILES = ["GSM2_barcodes.tsv.gz", "GSM2_features.tsv.gz", "GSM2_matrix.mtx.gz"]
IDAT_FILES = ["GSM3_201_R01C01_Grn.idat.gz", "GSM3_201_R01C01_Red.idat.gz"]


def _scrape(gse: str, files: list[str], title: str, tech: str = "") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Series": gse,
            "Platforms": "GPL1",
            "Platform_title": title,
            "Platform_technology": tech,
            "Supplementary file": files,
        }
    )


class DetectFileTechnologyTest(unittest.TestCase):
    def test_visium_beats_plain_10x(self):
        tech, matched = detect_file_technology(VISIUM_FILES)
        self.assertEqual(tech, "Visium")
        self.assertGreaterEqual(len(matched), 2)

    def test_signatures(self):
        self.assertEqual(detect_file_technology(XENIUM_FILES)[0], "Xenium")
        self.assertEqual(detect_file_technology(SCRNA_FILES)[0], "scRNA-seq")
        self.assertEqual(detect_file_technology(IDAT_FILES)[0], "DNA methylation")
        self.assertIsNone(detect_file_technology(["GSE1_counts.txt.gz"])[0])


class DeclaredTechnologyTest(unittest.TestCase):
    def test_mapping(self):
        self.assertEqual(declared_technology("Illumina NovaSeq 6000 (Homo sapiens)"), GENERIC)
        self.assertEqual(declared_technology("Xenium In Situ Analyzer"), "Xenium")
        self.assertEqual(
            declared_technology("Infinium MethylationEPIC", "oligonucleotide beads"),
            "DNA methylation",
        )


class FlagDubiousGsesTest(unittest.TestCase):
    def setUp(self):
        df = pd.concat(
            [
                _scrape("GSE277435", VISIUM_FILES, "Illumina NovaSeq 6000 (Homo sapiens)"),
                _scrape("GSE2", XENIUM_FILES, "Xenium In Situ Analyzer"),
                _scrape("GSE3", SCRNA_FILES, "Illumina NovaSeq 6000 (Homo sapiens)"),
                _scrape("GSE4", IDAT_FILES, "Infinium MethylationEPIC"),
                _scrape("GSE5", SCRNA_FILES, "Xenium In Situ Analyzer"),
            ],
            ignore_index=True,
        )
        self.result = flag_dubious_gses(df).set_index("Series")

    def test_visium_on_sequencer_is_generic_platform(self):
        row = self.result.loc["GSE277435"]
        self.assertTrue(row["Dubious"])
        self.assertEqual(row["File_signature_tech"], "Visium")
        self.assertEqual(row["Reason"], "generic platform")

    def test_matching_declarations_not_flagged(self):
        self.assertFalse(self.result.loc["GSE2", "Dubious"])
        self.assertFalse(self.result.loc["GSE4", "Dubious"])

    def test_scrna_on_sequencer_not_flagged(self):
        row = self.result.loc["GSE3"]
        self.assertEqual(row["File_signature_tech"], "scRNA-seq")
        self.assertEqual(row["Declared_tech"], GENERIC)
        self.assertFalse(row["Dubious"])

    def test_conflicting_platform(self):
        self.assertEqual(self.result.loc["GSE5", "Reason"], "conflicting platform")

    def test_empty_input(self):
        self.assertTrue(flag_dubious_gses(pd.DataFrame()).empty)


if __name__ == "__main__":
    unittest.main()
