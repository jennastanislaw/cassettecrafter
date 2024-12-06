import sys
import os
import pytest
import pandas as pd

# Tests
# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from mixed_base_rec_site_check import degen_codon_checker


# Mock Enzyme class for testing

# Updated Enzyme class
class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, oh_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = oh_length

    def __repr__(self):
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, oh_length={self.oh_length})"


# Test class with sample enzyme data
class TestReplaceEnzymeSite:
    @pytest.fixture
    def sample_enzyme_data(self):
        """Create a list of demo enzyme objects."""
        EcoRI = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        BamHI = Enzyme("BamHI", "GGATCC", "CCTAGG", 0, 6)
        FakeEnzyme = Enzyme("FakeEnzyme", "TGG", "TGG", 0, 3)

        return [EcoRI, BamHI, FakeEnzyme]

    # Integration Test 1: Basic Test Case (No Non-Canonical Bases)
    def test_no_non_canonical_bases(self, sample_enzyme_data):
        # Sample DataFrame with canonical bases
        test_df = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'DNA': ['ATGC', 'ATGC']
        })

        enzyme = sample_enzyme_data[0]  # Using EcoRI as the enzyme

        result_df = degen_codon_checker(test_df, enzyme)

        # Assert the output is identical to the input
        assert all(result_df['valid_mixed_bases'] == result_df['DNA'])

    # Integration Test 3: Test Case with Recognition Sites Overlap
    def test_recognition_site_overlap(self, sample_enzyme_data):
        df_with_overlap = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'DNA': ['ATGC', 'GAATNN']  # Sequence with a non-canonical base 'N'
        })

        enzyme = sample_enzyme_data[0]  # Using EcoRI as the enzyme

        test = enzyme.fwd_recognition_site

        result_df = degen_codon_checker(df_with_overlap, enzyme)

        # Assert that recognition site overlap is handled and added
        assert len(result_df) > len(df_with_overlap)
        assert not any(result_df['valid_mixed_bases'].str.contains('N'))  # Should contain expanded sequences with overlap

    # Integration Test 4: Test Case without Recognition Sites Overlap
    def test_recognition_site_no_overlap(self, sample_enzyme_data):
        df_no_overlap = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'DNA': ['ATGC', 'ATCN']  # Sequence with a non-canonical base 'N'
        })

        enzyme = sample_enzyme_data[2]  # Using FakeEnzyme (a non-canonical enzyme)

        result_df = degen_codon_checker(df_no_overlap, enzyme)

        # Assert that the original sequences are added as valid sequences
        assert all(
            result_df['valid_mixed_bases'] == result_df['DNA'])  # No expansion needed, original sequences are valid
