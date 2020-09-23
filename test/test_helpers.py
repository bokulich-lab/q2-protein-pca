import os
import numpy as np
import pandas as pd
from unittest import TestCase

from protein_pca.alignment import ProtAlign
from protein_pca.helpers import (
    split_string,
    calculate_covariance_matrix,
    euclidean_distance,
    make_label,
    positions_to_color,
    locate_positions,
)
from protein_pca.ranking import ProtRank


class TestHelpers(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.input_df = pd.DataFrame.from_dict(
            {
                "seq0": ["A", "B", "-", "C", "A"],
                "seq1": ["B", "-", "C", "C", "-"],
                "seq2": ["A", "B", "C", "D", "E"],
                "seq3": ["-", "A", "B", "C", "E"],
                "seq4": ["A", "B", "C", "D", "E"],
                "seq5": ["A", "B", "C", "E", "-"],
            },
            orient="index",
            columns=["pos1", "pos2", "pos3", "pos4", "pos5"],
        )
        cls.occurrence_table = pd.DataFrame.from_dict(
            {
                "-": [1, 1, 1, 0, 2],
                "A": [4, 1, 0, 0, 1],
                "B": [1, 4, 1, 0, 0],
                "C": [0, 0, 4, 3, 0],
                "D": [0, 0, 0, 2, 0],
                "E": [0, 0, 0, 1, 3],
            },
            orient="index",
            columns=["pos1", "pos2", "pos3", "pos4", "pos5"],
        )
        cls.cwd = os.path.dirname(os.path.realpath(__file__))
        cls.test_df = pd.read_csv(os.path.join(cls.cwd, "assets/test_data.csv"), index_col=0)
        cls.expected_ranks = pd.read_csv(
            os.path.join(cls.cwd, "assets/expected_ranks.csv"), index_col=0
        )
        cls.main_loc = os.path.join(cls.cwd, "test_main_loc")
        cls.results_loc = os.path.join(cls.cwd, "test_results")
        cls.ProtAlign = ProtAlign(
            main_location=cls.main_loc,
            results_location=cls.results_loc,
            input_fasta_file=os.path.join(cls.cwd, "assets/someseqs.fasta"),
            aln_slice_positions=os.path.join(cls.cwd, "alnpos.txt"),
            clustalo_location=None,
            clustalo_overwrite=True,
            seq_start=0,
            seq_end=-1,
        )
        cls.ProtRank = ProtRank(
            cls.main_loc, cls.results_loc, prot_align_obj=cls.ProtAlign
        )

    def test_split_string(self):
        input_str = "PROTEIN"
        expected_result = ["P", "R", "O", "T", "E", "I", "N"]
        actual_result = split_string(input_str)
        self.assertEqual(actual_result, expected_result)

    def test_get_occurrences(self):
        actual_result = self.ProtRank.get_occurrences(self.input_df)
        pd.testing.assert_frame_equal(actual_result, self.occurrence_table)

    def test_rank_alphabet(self):
        input_series = self.occurrence_table["pos5"]
        # TODO: double-check that the below result is correct
        expected_result = {"A": 1, "E": 2, "-": 0}
        actual_result = self.ProtRank.rank_alphabet(input_series)
        self.assertDictEqual(actual_result, expected_result)

        input_series = self.occurrence_table["pos4"]
        expected_result = {"C": 3, "D": 2, "E": 1, "-": 0}
        actual_result = self.ProtRank.rank_alphabet(input_series)
        self.assertDictEqual(actual_result, expected_result)

    def test_rank_all_columns(self):
        actual_result = self.ProtRank.rank_all_columns(self.test_df)
        pd.testing.assert_frame_equal(actual_result, self.expected_ranks)

    def test_variance_covariance(self):
        input_df = pd.DataFrame.from_dict(
            {
                1: [90, 60, 90],
                2: [90, 90, 30],
                3: [60, 60, 60],
                4: [60, 60, 90],
                5: [30, 30, 30],
            },
            orient="index",
            columns=["Math", "English", "Art"],
        )
        expected_result = np.array([[504, 360, 180], [360, 360, 0], [180, 0, 720]])
        actual_result = calculate_covariance_matrix(input_df)
        np.testing.assert_array_equal(actual_result, expected_result)

    def test_euclidean(self):
        expected_result = 2.236
        actual_result = euclidean_distance(1, 2)
        self.assertAlmostEqual(actual_result, expected_result, places=3)

        expected_result = 2.236
        actual_result = euclidean_distance(-1, 2)
        self.assertAlmostEqual(actual_result, expected_result, places=3)

        expected_result = 1
        actual_result = euclidean_distance(1, 0)
        self.assertAlmostEqual(actual_result, expected_result, places=3)

    def test_make_label(self):
        lbl = make_label(1)
        self.assertEqual(lbl, "0001")

        lbl = make_label(100)
        self.assertEqual(lbl, "0100")

        lbl = make_label(10000)
        self.assertEqual(lbl, "10000")

        lbl = make_label(10, trail_digits=6)
        self.assertEqual(lbl, "000010")

    def test_positions_to_color(self):
        position_list = [[2, 1.2, -2.4], [5, -2, 4]]
        expected_result = pd.DataFrame.from_dict(
            {0: [3, 1, "A", 1.2, 1.2], 1: [6, 4, "E", -2.0, 2.0]},
            orient="index",
            columns=["PosAln", "PosRaw", "AA", "PC score", "Abs PC score"],
        )
        actual_result = positions_to_color(position_list, self.input_df.iloc[3])
        pd.testing.assert_frame_equal(actual_result, expected_result)

    def test_locate_positions(self):
        x = [1, 2, 3, 4, 5]
        y = [0.5, 0.4, 0.3, 0.2, 0.1]
        expected_result = (
            [[0, 1, 0.5]],
            [[3, 4, 0.2], [4, 5, 0.1]],
            [[1, 2, 0.4], [2, 3, 0.3]],
        )
        actual_result = locate_positions(x, y, 0.4, 0.8)
        self.assertTupleEqual(actual_result, expected_result)
