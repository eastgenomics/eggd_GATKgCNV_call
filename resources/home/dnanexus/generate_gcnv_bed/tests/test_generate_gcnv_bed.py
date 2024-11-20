"""
Tests for generate_gcnv_bed.py
"""

from glob import glob
import os
from pathlib import Path
import shutil
from uuid import uuid4
import unittest
from unittest.mock import patch

import pandas as pd
import pytest

TEST_DATA_DIR = Path(__file__).absolute().parent.joinpath("test_data")

from ..generate_gcnv_bed import (
    calculate_mean_and_std_dev,
    read_single_copy_ratio_file,
    read_all_copy_ratio_files,
    write_run_level_bed_file,
    write_sample_bed_file,
)


class TestReadSingleCopyRatioFile(unittest.TestCase):
    sample_1_test_file = Path(TEST_DATA_DIR).joinpath(
        "sample_1_denoised_copy_ratios.tsv"
    )

    def test_sample_name_correctly_parsed_from_header(self):
        """
        Sample name is stored in `@RG` header line next to `SM:`, test
        that it is correctly read from the file
        """
        parsed_sample_name, _ = read_single_copy_ratio_file(
            copy_ratio_file=self.sample_1_test_file
        )

        self.assertEqual(parsed_sample_name, "sample_1")

    def test_contents_correctly_parsed_into_dataframe(self):
        expected_df = pd.DataFrame(
            {
                "chr": ["1", "1", "1", "1", "1"],
                "start": [7917094, 10566156, 17345176, 17345330, 17348949],
                "end": [7917213, 10566275, 17345329, 17345483, 17349114],
                "sample_1": [
                    1.9486853760805845,
                    2.0797460403307504,
                    1.9807613126697792,
                    2.0011478190013547,
                    1.9811431848626703,
                ],
            },
        )

        _, parsed_df = read_single_copy_ratio_file(
            copy_ratio_file=self.sample_1_test_file
        )

        self.assertTrue(parsed_df.equals(expected_df))

    def test_file_not_found_error_raised_when_invalid_file_path_given(self):
        with pytest.raises(
            FileNotFoundError,
            match="Copy ratio file not found: not_a_file.tsv",
        ):
            read_single_copy_ratio_file(copy_ratio_file="not_a_file.tsv")

    def test_value_error_raised_when_rg_line_not_present_in_header(self):
        with pytest.raises(
            ValueError,
            match=f"'@RG' header not found in {Path(__file__).absolute()}",
        ):
            # test we catch missing @RG header by passing in the test script
            # instead of value copy ratio file
            read_single_copy_ratio_file(
                copy_ratio_file=Path(__file__).absolute()
            )

    def test_value_error_raised_when_empty_file_provided_to_read_from(self):
        empty_test_file = f"{uuid4().hex}.empty.txt"

        with open(empty_test_file, encoding="utf-8", mode="w") as fh:
            # generate output file with only header and no data
            fh.write("@RG	ID:GATKCopyNumber	SM:sample_1")

        with self.subTest():
            with pytest.raises(
                ValueError, match=f"No data found in file: {empty_test_file}"
            ):
                read_single_copy_ratio_file(copy_ratio_file=empty_test_file)

        os.remove(empty_test_file)

    @patch("generate_gcnv_bed.generate_gcnv_bed.pd.read_csv")
    def test_value_error_raised_when_error_raised_from_pd_read_csv(
        self, mock_read
    ):
        """
        Ensure we catch any error that may be raised from pd.read_csv
        """
        mock_read.side_effect = RuntimeError()

        with pytest.raises(
            ValueError, match=f"Error reading file {self.sample_1_test_file}"
        ):
            read_single_copy_ratio_file(
                copy_ratio_file=self.sample_1_test_file
            )


class TestReadAllCopyRatioFiles(unittest.TestCase):
    def test_reading_all_files_to_dataframe_correct(self):
        copy_ratio_files = glob(f"{TEST_DATA_DIR}/*copy_ratios.tsv")

        copy_ratio_df = read_all_copy_ratio_files(
            copy_ratio_files=copy_ratio_files
        )

        expected_df = pd.DataFrame(
            {
                "chr": ["1", "1", "1", "1", "1"],
                "start": [7917094, 10566156, 17345176, 17345330, 17348949],
                "end": [7917213, 10566275, 17345329, 17345483, 17349114],
                "sample_1": [
                    1.9486853760805845,
                    2.0797460403307504,
                    1.9807613126697792,
                    2.0011478190013547,
                    1.9811431848626703,
                ],
                "sample_2": [
                    1.6825420819612993,
                    1.8989351264973224,
                    1.949430403226796,
                    1.93894585147577,
                    2.015844320949207,
                ],
                "sample_3": [
                    2.0369879930267016,
                    2.0132574495419986,
                    1.9679869323497048,
                    2.139656851102916,
                    1.9622012727554363,
                ],
            },
        )

        self.assertTrue(copy_ratio_df.equals(expected_df))

    def test_value_error_raised_when_file_with_different_intervals_parsed(
        self,
    ):
        """
        Test we correctly raise a ValueError when not all files have the
        same intervals by making an additional test copy ratio file
        and dropping some rows
        """
        diff_intervals_tsv = (
            Path(__file__)
            .absolute()
            .parent.joinpath("diff_intervals_copy_ratios.tsv")
        )

        with open(
            Path(TEST_DATA_DIR).joinpath("sample_1_denoised_copy_ratios.tsv"),
            encoding="utf-8",
            mode="r",
        ) as fh:
            sample_1_contents = fh.read().splitlines()

        with open(
            diff_intervals_tsv,
            encoding="utf-8",
            mode="w",
        ) as fh:
            fh.write("\n".join(sample_1_contents[:-2]))

        copy_ratio_files = glob(f"{TEST_DATA_DIR}/*copy_ratios.tsv")
        copy_ratio_files.append(diff_intervals_tsv)

        with self.subTest():
            with pytest.raises(
                ValueError,
                match="Copy ratio file for sample_1",
            ):
                read_all_copy_ratio_files(copy_ratio_files=copy_ratio_files)

        os.remove(diff_intervals_tsv)
