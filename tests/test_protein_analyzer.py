

import unittest
from pathlib import Path
import sys
from collections import Counter
from types import SimpleNamespace

# Import the actual Protein Analyzer
PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))

import protein_analyzer as pa


def get_args():
    """Create default arguments for protein analysis."""
    return SimpleNamespace(
        ph=7.4,
        force=False,
        hydropathy=False,
        window=19,
        hydro_threshold=1.6
    )


class TestProteinAnalyzer(unittest.TestCase):

    def test_protein_length(self):
        result = pa.analyze_record(
            "TestProtein", "ACDEFGHIK", get_args()
        )
        self.assertEqual(result["length"], 9)

    def test_amino_acid_composition(self):
        counts = Counter("AAACCG")

        self.assertEqual(counts["A"], 3)
        self.assertEqual(counts["C"], 2)
        self.assertEqual(counts["G"], 1)

    def test_molecular_weight(self):
        counts = Counter("ACDEFGHIK")

        minimum, maximum = pa.molecular_weight_range(counts)

        self.assertAlmostEqual(minimum, 1019.14004, places=3)
        self.assertEqual(minimum, maximum)

    def test_isoelectric_point(self):
        sequence = "ACDEFGHIK"
        counts = Counter(sequence)

        pi = pa.isoelectric_point(
            counts, sequence[0], sequence[-1]
        )

        self.assertGreaterEqual(pi, 0)
        self.assertLessEqual(pi, 14)

        charge_at_pi = pa.net_charge(
            counts, sequence[0], sequence[-1], pi
        )

        self.assertAlmostEqual(charge_at_pi, 0, places=3)

    def test_gravy_calculation(self):
        score, skipped = pa.gravy("ACDEFGHIK")

        self.assertAlmostEqual(score, -0.322222, places=5)
        self.assertEqual(skipped, 0)

    def test_ambiguous_amino_acid(self):
        counts = Counter("ACDX")

        minimum, maximum = pa.molecular_weight_range(counts)

        self.assertLess(minimum, maximum)

    def test_invalid_amino_acid_removal(self):
        result = pa.analyze_record(
            "InvalidTest", "ACD123", get_args()
        )

        self.assertEqual(result["length"], 3)
        self.assertEqual(result["invalid"], 3)
        self.assertFalse(result["skipped"])

    def test_empty_sequence(self):
        result = pa.analyze_record(
            "EmptyTest", "", get_args()
        )

        self.assertTrue(result["skipped"])
        self.assertEqual(result["length"], 0)

    def test_nucleotide_sequence_detection(self):
        dna = "ATGCGTACGTAGCTAGCTAGCTAG"

        result = pa.analyze_record(
            "DNA_Test", dna, get_args()
        )

        self.assertTrue(result["skipped"])

    def test_instability_index(self):
        score = pa.instability_index("ACDEFGHIK")

        self.assertIsInstance(score, float)


if __name__ == "__main__":
    unittest.main()
