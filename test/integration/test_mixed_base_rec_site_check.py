"""
This script defines tests for the `degen_codon_checker` function from the `mixed_base_rec_site_check` module.
The purpose of the function is to check and handle degenerate codons (IUPAC codes) within DNA sequences, particularly
in relation to enzyme recognition sites. The script contains the following components:

1. **Enzyme Class**:
   - Represents an enzyme with properties such as its name, forward and reverse recognition sites,
     spacer length, and overhang length. This class is used to provide sample data for the tests.

2. **TestReplaceEnzymeSite Class**:
   - Contains several test methods to validate the functionality of the `degen_codon_checker` function.
   - Each test method prepares a DataFrame of DNA sequences, applies the `degen_codon_checker` function,
     and verifies the output according to different scenarios, such as handling IUPAC codes, recognition site overlap,
     and no overlap.
"""
import sys
import os
import pytest
import pandas as pd

# Add the src directory to the Python path to allow importing the 'degen_codon_checker' function.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from mixed_base_rec_site_check import degen_codon_checker


class Enzyme:
    """
    Represents an enzyme with its name, forward and reverse recognition sites,
    spacer length, and overhang length. This class is used for enzyme data in the tests.

    Attributes:
        name (str): The name of the enzyme.
        fwd_recognition_site (str): The forward recognition site of the enzyme.
        rev_recognition_site (str): The reverse recognition site of the enzyme.
        spacer_length (int): The length of the spacer sequence.
        OH_length (int): The length of the overhang sequence.
    """

    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, oh_length):
        """
        Initializes the Enzyme object with the provided parameters.

        Args:
            name (str): The name of the enzyme.
            fwd_recognition_site (str): The forward recognition site of the enzyme.
            rev_recognition_site (str): The reverse recognition site of the enzyme.
            spacer_length (int): The length of the spacer sequence.
            oh_length (int): The length of the overhang sequence.
        """
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = oh_length

    def __repr__(self):
        """Returns a string representation of the Enzyme object."""
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, oh_length={self.OH_length})"


class TestReplaceEnzymeSite:
    """
    A set of unit tests for the `degen_codon_checker` function, which checks and handles degenerate codons
    (IUPAC codes) in DNA sequences relative to enzyme recognition sites.

    The tests cover several scenarios such as handling IUPAC codes, managing overlapping recognition sites,
    and checking sequences without the need for modification.
    """

    @pytest.fixture
    def sample_enzyme_data(self):
        """
        Creates a list of demo enzyme objects for use in the tests. This fixture provides sample data
        that simulates real enzyme recognition sites.

        Returns:
            list: A list of Enzyme objects, including EcoRI, BamHI, and FakeEnzyme.
        """
        EcoRI = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        BamHI = Enzyme("BamHI", "GGATCC", "CCTAGG", 0, 6)
        FakeEnzyme = Enzyme("FakeEnzyme", "TGG", "TGG", 0, 3)

        return [EcoRI, BamHI, FakeEnzyme]

    def test_no_IUPAC_codes(self, sample_enzyme_data):
        """
        Test when there are no IUPAC codes in the DNA sequences.
        The test ensures that the input DNA remains unchanged as no modification is required.

        Args:
            sample_enzyme_data (list): A list of Enzyme objects used for testing.

        Asserts:
            The 'valid_mixed_bases' column in the output DataFrame should match the original 'DNA' sequences.
        """
        # Sample DataFrame with valid DNA sequences (no IUPAC codes)
        test_df = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'DNA': ['ATGC', 'ATGC']
        })

        enzyme = sample_enzyme_data[0]  # Using EcoRI as the enzyme

        result_df = degen_codon_checker(test_df, enzyme)

        # Assert that the output is identical to the input
        assert all(result_df['valid_mixed_bases'] == result_df['DNA'])

    def test_recognition_site_overlap(self, sample_enzyme_data):
        """
        Test when the DNA sequences contain an IUPAC code ('N') that overlaps with the enzyme's recognition site.
        The test checks that the function correctly handles the overlap and expands the sequences.

        Args:
            sample_enzyme_data (list): A list of Enzyme objects used for testing.

        Asserts:
            The result should have no 'N' in the 'valid_mixed_bases' and the length of the output should be greater
            than the input (indicating that sequences were expanded due to the overlap).
        """
        df_with_overlap = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'DNA': ['ATGC', 'GAATNN']  # Sequence with IUPAC code 'N'
        })

        enzyme = sample_enzyme_data[0]  # Using EcoRI as the enzyme

        result_df = degen_codon_checker(df_with_overlap, enzyme)

        # Assert that recognition site overlap is handled and added
        assert len(result_df) > len(df_with_overlap)
        assert not any(
            result_df['valid_mixed_bases'].str.contains('N'))  # Should contain expanded sequences with overlap

    def test_recognition_site_no_overlap(self, sample_enzyme_data):
        """
        Test when the DNA sequences contain an IUPAC code ('N') but do not overlap with the enzyme's recognition site.
        The test ensures that no modification is made if no overlap occurs.

        Args:
            sample_enzyme_data (list): A list of Enzyme objects used for testing.

        Asserts:
            The original sequences should remain unchanged in the output DataFrame.
        """
        df_no_overlap = pd.DataFrame({
            'name': ['seq1', 'seq2'],
            'DNA': ['ATGC', 'ATCN']  # Sequence with IUPAC code 'N'
        })

        enzyme = sample_enzyme_data[2]  # Using FakeEnzyme (a non-canonical enzyme)

        result_df = degen_codon_checker(df_no_overlap, enzyme)

        # Assert that the original sequences are added as valid sequences
        assert all(
            result_df['valid_mixed_bases'] == result_df['DNA'])  # No expansion needed, original sequences are valid
