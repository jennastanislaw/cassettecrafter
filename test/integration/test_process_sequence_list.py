"""
This script contains a series of test cases to validate the behavior of the 'process_dna_sequences'
function, which processes DNA sequences based on the recognition sites of a given restriction enzyme.
The test suite checks various scenarios to ensure that the function behaves as expected under different
conditions.

Test cases include:
1. **test_basic_case**: Verifies that sequences contain exactly one instance of the enzyme's forward
   and reverse recognition sites, and that sequences meet the minimum oligo size requirement.
2. **test_random_oligos**: Ensures that the function extends sequences to meet the minimum oligo size
   when necessary.
3. **test_no_cassette_columns**: Checks the behavior when no cassette columns are present in the
   dataframe, ensuring the dataframe remains unchanged.
4. **test_empty_dataframe**: Tests the function's behavior when an empty dataframe is provided, ensuring
   that the result is also an empty dataframe.

Fixtures:
- **enzyme**: Provides a test enzyme object with specific recognition sites, spacer length, and
  overhang length.
- **sample_dataframe**: Provides a sample dataframe with cassette DNA sequences for testing.

The tests aim to ensure that the 'process_dna_sequences' function handles various edge cases and
processes DNA sequences as expected, taking into account enzyme recognition sites and sequence length
constraints.
"""

import pytest
import pandas as pd
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from process_sequence_list import process_dna_sequences


class Enzyme:
    """
    A class representing a restriction enzyme with its recognition sites and other relevant properties.
    """

    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        """
        Initializes an Enzyme object with the given properties.

        Parameters:
        - name (str): The name of the enzyme.
        - fwd_recognition_site (str): The forward recognition site sequence of the enzyme.
        - rev_recognition_site (str): The reverse recognition site sequence of the enzyme.
        - spacer_length (int): The length of the spacer sequence for the enzyme.
        - OH_length (int): The length of the overhang generated after digestion.
        """
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = OH_length

    def __repr__(self):
        """
        Returns a string representation of the enzyme object.

        Returns:
        - str: A string describing the enzyme and its properties.
        """
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, OH_length={self.OH_length})"


class TestProcessDnaSequences:
    """
    Test suite for the process_dna_sequences function, which processes DNA sequences
    by applying an enzyme's recognition sites and ensuring that they meet certain criteria.
    """

    @pytest.fixture
    def enzyme(self):
        """
        Fixture to provide a test enzyme object.

        This enzyme has a forward recognition site "GAATTC", reverse recognition site "CTTAAG",
        a spacer length of 0, and an overhang length of 6.

        Returns:
        - Enzyme: A test enzyme object.
        """
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", spacer_length=0, OH_length=6)

    @pytest.fixture
    def sample_dataframe(self):
        """
        Fixture to provide a sample dataframe for testing.

        The dataframe contains two cassette columns ("Cassette_1" and "Cassette_2") with sample DNA sequences.

        Returns:
        - pd.DataFrame: A dataframe with sample DNA sequences.
        """
        data = {
            "Cassette_1": ["ATGCGT", "CGTAGC", "GCTAGC"],
            "Cassette_2": ["GATCGA", "TACGTA", "CGATCG"],
        }
        return pd.DataFrame(data)

    def test_basic_case(self, sample_dataframe, enzyme):
        """
        Test the process_dna_sequences function with basic input data.

        This test verifies that the sequences contain exactly one instance of the enzyme's forward and reverse
        recognition sites and that all sequences meet the minimum oligo size.
        """
        min_oligo_size = 10  # Set the minimum oligo size for this test

        # Process the DNA sequences using the enzyme and the minimum oligo size
        result_df = process_dna_sequences(sample_dataframe, enzyme, min_oligo_size)

        # Verify that each sequence contains exactly one instance of the forward and reverse recognition sites
        for col in ["Cassette_1", "Cassette_2"]:
            assert all(seq.count(enzyme.fwd_recognition_site) == 1 for seq in result_df[col]), \
                f"Some sequences in {col} do not contain exactly one instance of the forward recognition site."
            assert all(seq.count(enzyme.rev_recognition_site) == 1 for seq in result_df[col]), \
                f"Some sequences in {col} do not contain exactly one instance of the reverse recognition site."

        # Verify that all sequences meet the minimum oligo size
        assert all(result_df["Cassette_1"].str.len() >= min_oligo_size), \
            "Some sequences in Cassette_1 do not meet the minimum oligo size."
        assert all(result_df["Cassette_2"].str.len() >= min_oligo_size), \
            "Some sequences in Cassette_2 do not meet the minimum oligo size."

    def test_random_oligos(self, sample_dataframe, enzyme):
        """
        Test the process_dna_sequences function when random oligos need to be generated.

        This test checks if the sequences are extended to meet the minimum oligo size when required.
        """
        min_oligo_size = 20  # Set a higher threshold to force random oligo generation

        # Process the DNA sequences and check if the sequences meet the minimum size
        result_df = process_dna_sequences(sample_dataframe, enzyme, min_oligo_size)

        # Verify that sequences have been extended to meet the minimum size
        for col in ["Cassette_1", "Cassette_2"]:
            assert all(result_df[col].str.len() >= min_oligo_size)

    def test_no_cassette_columns(self, enzyme):
        """
        Test the behavior when no cassette columns are present in the dataframe.

        This test verifies that the dataframe remains unchanged if no valid cassette columns are found.
        """
        df = pd.DataFrame({"NonCassette_1": ["ATGCGT"], "NonCassette_2": ["GATCGA"]})
        min_oligo_size = 10  # Set the minimum oligo size

        # Process the DNA sequences and ensure the dataframe remains unchanged
        result_df = process_dna_sequences(df, enzyme, min_oligo_size)

        # Verify that the dataframe is unchanged
        pd.testing.assert_frame_equal(df, result_df)

    def test_empty_dataframe(self, enzyme):
        """
        Test the behavior when an empty dataframe is provided.

        This test checks if the function handles empty dataframes without errors.
        """
        df = pd.DataFrame()  # Create an empty dataframe
        min_oligo_size = 10  # Set the minimum oligo size

        # Process the empty dataframe and verify that it remains empty
        result_df = process_dna_sequences(df, enzyme, min_oligo_size)

        # Verify that the dataframe remains empty
        assert result_df.empty
