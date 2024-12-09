"""
This script defines a set of tests for the `replace_enzyme_site` function, which is responsible for replacing
specific recognition sites in a DNA sequence based on enzyme properties. The script includes the following components:

1. **Enzyme Class**:
   - A class representing an enzyme with attributes such as its name, forward and reverse recognition sites,
     spacer length, and overhang length.

2. **TestReplaceEnzymeSite Class**:
   - A test suite that verifies the behavior of the `replace_enzyme_site` function. The test suite includes several
     test methods to cover different scenarios:
     - **Successful site replacement**: Verifies that a recognition site is successfully replaced with a new codon.
     - **No replacement possible**: Tests the case where no suitable replacement can be found for a given recognition site.
     - **Multiple sites replacement**: Ensures that multiple recognition sites in a DNA sequence can be replaced correctly.
     - **Empty DNA sequence**: Verifies the behavior when an empty DNA sequence is passed to the function.

The tests use the `pytest` framework, with a fixture providing demo enzyme data and methods ensuring that
the `replace_enzyme_site` function behaves as expected in various scenarios.

The script is primarily focused on testing the enzyme site replacement functionality, ensuring that the
logic works as intended under different conditions and edge cases.
"""

import sys
import os
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
from enzyme_site_replacement import replace_enzyme_site


class Enzyme:
    """
    Represents an enzyme with recognition sites and overhang properties.

    Attributes:
        name (str): The name of the enzyme.
        fwd_recognition_site (str): The forward recognition sequence of the enzyme.
        rev_recognition_site (str): The reverse recognition sequence of the enzyme.
        spacer_length (int): The length of the spacer (non-recognition) sequence.
        oh_length (int): The length of the overhang sequence.

    Methods:
        __repr__(): Returns a string representation of the enzyme object.
    """

    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, oh_length):
        """
        Initializes the Enzyme object with the provided attributes.

        Args:
            name (str): The name of the enzyme.
            fwd_recognition_site (str): The forward recognition site.
            rev_recognition_site (str): The reverse recognition site.
            spacer_length (int): The length of the spacer sequence.
            oh_length (int): The length of the overhang sequence.
        """
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.oh_length = oh_length

    def __repr__(self):
        """
        Returns a string representation of the enzyme object.

        Returns:
            str: The string representation of the enzyme object with its attributes.
        """
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, oh_length={self.oh_length})"


class TestReplaceEnzymeSite:
    """
    A test suite for the `replace_enzyme_site` function, which replaces recognition sites in a given DNA sequence
    with suitable alternatives based on enzyme properties.

    The tests cover various cases including:
    - Successful site replacement.
    - Situations where no replacement is possible.
    - Replacement for multiple recognition sites.
    - Handling of empty DNA sequences.
    """

    @pytest.fixture
    def sample_enzyme_data(self):
        """
        Creates a list of demo enzyme objects to be used in the tests.

        Returns:
            list: A list of Enzyme objects for testing purposes.
        """
        EcoRI = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        BamHI = Enzyme("BamHI", "GGATCC", "CCTAGG", 0, 4)
        FakeEnzyme = Enzyme("FakeEnzyme", "TGG", "TGG", 0, 4)

        return [EcoRI, BamHI, FakeEnzyme]

    def test_successful_replacement(self, sample_enzyme_data):
        """
        Tests the successful replacement of a recognition site in a DNA sequence.

        Specifically, this test checks if the EcoRI recognition site ("GAATTC") is replaced
        with a new codon while maintaining the sequence functionality.
        """
        # Extract EcoRI enzyme from sample data
        EcoRI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "EcoRI")
        DNA = "AGCTGAATTCAGCTCTTAAGCTAGATG"  # DNA sequence containing EcoRI recognition site
        new_DNA = replace_enzyme_site(EcoRI, DNA)

        # Assert that the EcoRI recognition site is not present in the new DNA
        assert "GAATTC" not in new_DNA
        assert "GAA" in new_DNA or "TTC" in new_DNA  # Check for synonymous codons after replacement

    def test_no_replacement_possible(self, sample_enzyme_data):
        """
        Tests the scenario where no replacement is possible due to the lack of suitable
        replacement options for a given recognition site.

        Specifically, it checks if the function raises a ValueError when trying to replace
        a non-functional enzyme's recognition site.
        """
        # Extract FakeEnzyme from sample data
        FakeEnzyme = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "FakeEnzyme")
        DNA = "TGGTGAGATCCAGCTCTTAAGCTAGATG"  # DNA sequence with BamHI recognition site
        key = "0"  # Expected failure because FakeEnzyme has no suitable recognition site

        # Assert that a ValueError is raised indicating no replacement was found
        with pytest.raises(ValueError, match=f"No suitable replacement found for recognition site at index {key}."):
            replace_enzyme_site(FakeEnzyme, DNA)

    def test_multiple_sites_replacement(self, sample_enzyme_data):
        """
        Tests the replacement of multiple recognition sites in a DNA sequence.

        This test covers a case where both EcoRI and BamHI recognition sites exist in the
        sequence, and ensures both sites are successfully replaced.
        """
        # Extract EcoRI enzyme from sample data
        EcoRI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "EcoRI")
        DNA = "AGCTGAATTCGGATCCGCTCTTAAGCTAGATG"  # DNA with both EcoRI and BamHI recognition sites
        new_DNA = replace_enzyme_site(EcoRI, DNA)

        # Assert that EcoRI site is replaced
        assert "GAATTC" not in new_DNA

        # Extract BamHI enzyme from sample data
        BamHI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "BamHI")
        new_DNA = replace_enzyme_site(BamHI, new_DNA)

        # Assert that BamHI site is replaced
        assert "GGATCC" not in new_DNA

    def test_empty_dna_sequence(self, sample_enzyme_data):
        """
        Tests the behavior of the `replace_enzyme_site` function when given an empty DNA sequence.

        It checks that the function correctly returns an empty string when no DNA sequence is provided.
        """
        # Extract EcoRI enzyme from sample data
        EcoRI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "EcoRI")
        DNA = ""  # Empty DNA sequence
        new_DNA = replace_enzyme_site(EcoRI, DNA)

        # Assert that the new DNA is still an empty string
        assert new_DNA == ""
