import tempfile
import unittest
from pathlib import Path

from predict import BLAST_COLS, read_blast


class BlastSchemaTest(unittest.TestCase):
    def test_reads_legacy_and_esm_rows(self):
        row = ["0"] * len(BLAST_COLS)
        row[0:5] = ["chr1", "10", "20", "read_ORF1", "+"]
        row[-1] = "UNKNOWN"

        with tempfile.TemporaryDirectory() as directory:
            legacy = Path(directory) / "legacy.tsv"
            esm = Path(directory) / "esm.tsv"
            legacy.write_text("\t".join(row) + "\n", encoding="utf-8")
            esm.write_text("\t".join(row[:-1] + ["73.5", row[-1]]) + "\n", encoding="utf-8")

            self.assertNotIn("esm_plddt", read_blast(legacy).columns)
            self.assertEqual(read_blast(esm).loc[0, "esm_plddt"], 73.5)


if __name__ == "__main__":
    unittest.main()
