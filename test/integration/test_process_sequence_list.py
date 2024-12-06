import pytest
import pandas as pd
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from process_sequence_list import process_dna_sequences

class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = OH_length

    def __repr__(self):
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, OH_length={self.OH_length})"

class TestProcessDnaSequences:
    """Test suite for the process_dna_sequences function."""

    @pytest.fixture
    def enzyme(self):
        """Fixture to provide a test enzyme."""
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", spacer_length=0, OH_length=6)

    @pytest.fixture
    def sample_dataframe(self):
        """Fixture to provide a sample dataframe for testing."""
        data = {
            "Cassette_1": ["ATGCGT", "CGTAGC", "GCTAGC"],
            "Cassette_2": ["GATCGA", "TACGTA", "CGATCG"],
        }
        return pd.DataFrame(data)

    def test_basic_case(self, sample_dataframe, enzyme):
        """Test the function with basic input."""
        min_oligo_size = 10

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
        """Test the function when random oligos need to be generated."""
        min_oligo_size = 20  # Higher threshold to force random oligo generation

        result_df = process_dna_sequences(sample_dataframe, enzyme, min_oligo_size)

        # Verify that sequences have been extended to meet the minimum size
        for col in ["Cassette_1", "Cassette_2"]:
            assert all(result_df[col].str.len() >= min_oligo_size)

    def test_no_cassette_columns(self, enzyme):
        """Test behavior when no cassette columns are present."""
        df = pd.DataFrame({"NonCassette_1": ["ATGCGT"], "NonCassette_2": ["GATCGA"]})
        min_oligo_size = 10

        result_df = process_dna_sequences(df, enzyme, min_oligo_size)

        # Verify that the dataframe is unchanged
        pd.testing.assert_frame_equal(df, result_df)

    def test_empty_dataframe(self, enzyme):
        """Test behavior with an empty dataframe."""
        df = pd.DataFrame()
        min_oligo_size = 10

        result_df = process_dna_sequences(df, enzyme, min_oligo_size)

        # Verify that the dataframe remains empty
        assert result_df.empty

