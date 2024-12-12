"""
This script contains tests for functions that handle enzyme recognition site
detection and sequence validation for DNA sequence manipulation. The tests
cover the functionality of recognizing enzyme recognition sites, handling
IUPAC codes, and validating sequences according to specified rules.

Dependencies:
- `sys`: For modifying the Python path.
- `os`: For handling file paths and directories.
- `enzyme_site_replacement_utils`: Contains functions for enzyme recognition
  and sequence validation.
- `pytest`: For running and asserting test cases.
"""

import sys
import os
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from mixed_base_rec_site_check_utils import (
    expand_dna_sequence,
    find_IUPAC_codes,
    check_IUPAC_code_in_rec_sites,
    check_recognition_sites_in_expanded_sequences,
    append_valid_sequences
)
class Enzyme:
    """
    A class representing a restriction enzyme with its properties, including
    the recognition sites (forward and reverse), spacer length, and overhang length.

    Attributes:
        name (str): The name of the enzyme.
        fwd_recognition_site (str): The forward recognition site sequence.
        rev_recognition_site (str): The reverse recognition site sequence.
        spacer_length (int): The length of the spacer sequence between the recognition site and the cutting site.
        OH_length (int): The overhang length, which is the length of the DNA fragment left after the enzyme cuts the DNA.

    Methods:
        __repr__(): Returns a string representation of the Enzyme object, useful for debugging and inspecting enzyme attributes.
    """

    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, oh_length):
        """
        Initializes a new instance of the Enzyme class with the given attributes.

        Parameters:
            name (str): The name of the enzyme.
            fwd_recognition_site (str): The forward recognition site sequence.
            rev_recognition_site (str): The reverse recognition site sequence.
            spacer_length (int): The length of the spacer sequence between the recognition site and the cutting site.
            oh_length (int): The overhang length of the DNA fragment after the enzyme cuts the DNA.
        """
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = oh_length

    def __repr__(self):
        """
        Returns a string representation of the Enzyme object.

        The returned string contains the enzyme name, its forward and reverse recognition sites,
        spacer length, and overhang length in a readable format. This is useful for debugging
        and inspecting enzyme attributes.

        Returns:
            str: A string representation of the Enzyme object.
        """
        return (f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, "
                f"rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, "
                f"oh_length={self.OH_length})")


class TestExpandDnaSequence:
    """
    This class contains tests for the expand_dna_sequence function, which
    handles the expansion of DNA sequences containing IUPAC codes. IUPAC codes
    represent ambiguity in DNA bases and allow multiple possible base pairings.
    The tests ensure that the function correctly expands sequences containing
    valid IUPAC codes, handles invalid bases, and processes edge cases.
    """

    def test_single_base(self):
        """
        Test the expansion of single bases, including valid IUPAC codes.
        - 'A' should return ['A'] as it represents adenine.
        - 'N' should return ['A', 'C', 'G', 'T'] as 'N' can be any of the four bases.
        - 'R' should return ['A', 'G'] as 'R' represents purine, which includes adenine and guanine.
        """
        assert expand_dna_sequence("A") == ["A"]  # A valid single base.
        assert expand_dna_sequence("N") == ["A", "C", "G", "T"]  # N can be any base.
        assert expand_dna_sequence("R") == ["A", "G"]  # R represents purines (A or G).

    def test_multiple_bases(self):
        """
        Test sequences containing multiple valid IUPAC codes.
        - 'AC' should return ['AC'] as it contains no ambiguity.
        - 'AN' should return all combinations of 'A' and 'N' (A, C, G, T).
        - 'NR' should return all possible combinations between 'N' and 'R'.
        """
        assert expand_dna_sequence("AC") == ["AC"]  # No ambiguity, should return as is.
        assert sorted(expand_dna_sequence("AN")) == ["AA", "AC", "AG", "AT"]  # All combinations of A and N.
        # All combinations of N and R.
        assert sorted(expand_dna_sequence("NR")) == ["AA", "AG", "CA", "CG", "GA", "GG", "TA", "TG"]

    def test_invalid_base(self):
        """
        Test sequences containing invalid bases (e.g., bases that are not in the IUPAC code).
        - A ValueError should be raised if an invalid base is encountered.
        """
        with pytest.raises(ValueError, match="Invalid base 'X' in DNA sequence."):
            expand_dna_sequence("AX")  # 'X' is not a valid base.

    def test_empty_sequence(self):
        """
        Test the case where the input DNA sequence is empty.
        - An empty input should return a list containing an empty string.
        """
        assert expand_dna_sequence("") == [""]  # Empty input should return [''].

    def test_complex_sequence(self):
        """
        Test a complex DNA sequence with multiple IUPAC codes.
        - The function should correctly expand a sequence containing multiple IUPAC codes.
        - In this case, 'RYN' should expand to a list of all combinations of 'R', 'Y', and 'N'.
        """
        dna_sequence = "RYN"
        expected = [
            "ACA", "ACC", "ACG", "ACT",
            "ATA", "ATC", "ATG", "ATT",
            "GCA", "GCC", "GCG", "GCT",
            "GTA", "GTC", "GTG", "GTT",
        ]
        # Check if the expanded sequence matches the expected output.
        assert sorted(expand_dna_sequence(dna_sequence)) == sorted(expected)


class TestFindIUPACcodes:
    """
    This class contains tests for the find_IUPAC_codes function, which identifies
    IUPAC codes (ambiguous nucleotide bases) in a DNA sequence. The function
    returns the positions and corresponding IUPAC codes found in the sequence.
    The tests cover cases with no IUPAC codes, all IUPAC codes, mixed sequences,
    and sequences with special characters or lowercase bases.
    """

    def test_no_IUPAC_codes(self):
        """
        Test a sequence with no IUPAC codes (standard DNA bases only).
        - The function should return an empty list if no IUPAC codes are found.
        """
        sequence = "ACGTACGT"  # No IUPAC codes in this sequence.
        expected = []  # No IUPAC codes found.
        assert find_IUPAC_codes(sequence) == expected  # Assert that the result matches the expected output.

    def test_all_IUPAC_codes(self):
        """
        Test a sequence consisting solely of IUPAC codes.
        - The function should return a list of tuples, each containing the index and corresponding IUPAC code.
        """
        sequence = "RYWSKMBDHVN"  # A sequence with only IUPAC codes.
        expected = [
            (0, 'R'), (1, 'Y'), (2, 'W'), (3, 'S'), (4, 'K'),
            (5, 'M'), (6, 'B'), (7, 'D'), (8, 'H'), (9, 'V'), (10, 'N')
        ]  # Expected positions and IUPAC codes.
        assert find_IUPAC_codes(sequence) == expected  # Assert the function returns the correct list of tuples.

    def test_mixed_sequence(self):
        """
        Test a sequence containing both standard DNA bases and IUPAC codes.
        - The function should identify only the positions of IUPAC codes.
        """
        sequence = "ACGTNACGTRY"  # Sequence with some IUPAC codes ('N', 'R', 'Y').
        expected = [(4, 'N'), (9, 'R'), (10, 'Y')]  # Expected positions and IUPAC codes.
        assert find_IUPAC_codes(sequence) == expected  # Assert the function correctly identifies IUPAC codes.

    def test_empty_sequence(self):
        """
        Test an empty DNA sequence.
        - The function should return an empty list for an empty input.
        """
        sequence = ""  # Empty sequence.
        expected = []  # No IUPAC codes in an empty sequence.
        assert find_IUPAC_codes(sequence) == expected  # Assert that the result is an empty list.

    def test_lowercase_bases(self):
        """
        Test a sequence with lowercase bases.
        - The function should be case-insensitive and identify IUPAC codes in lowercase as well.
        """
        sequence = "acgtnRY"  # Sequence with lowercase bases and IUPAC codes.
        expected = [(4, 'n'), (5, 'R'), (6, 'Y')]  # Expected positions and IUPAC codes, case-insensitive.
        assert find_IUPAC_codes(sequence) == expected  # Assert the function identifies IUPAC codes regardless of case.

    def test_special_characters(self):
        """
        Test a sequence containing special characters that are not IUPAC codes.
        - The function should identify and return positions of non-standard characters in the sequence.
        """
        sequence = "ACGT@!NRY"  # Sequence with special characters '@' and '!', and IUPAC codes.
        expected = [(4, '@'), (5, '!'), (6, 'N'), (7, 'R'), (8, 'Y')]  # Expected positions of special characters and IUPAC codes.
        assert find_IUPAC_codes(sequence) == expected  # Assert that the function correctly identifies special characters and IUPAC codes.


class TestCheckIUPACCodeInRecSites:
    """
    This class contains tests for the check_IUPAC_code_in_rec_sites function, which checks
    whether IUPAC codes (ambiguous nucleotide bases) are located within specified recognition
    sites (forward and reverse) based on the enzyme overhang length. The function returns
    True if any IUPAC code overlaps with a recognition site, and False otherwise.
    The tests cover a range of cases such as no IUPAC codes, IUPAC codes inside and outside
    recognition sites, and edge cases where IUPAC codes are at the boundary of recognition sites.
    """

    def test_no_IUPAC_codes(self):
        """
        Test case where no IUPAC codes are present in the sequence.
        - The function should return False, as no IUPAC codes can overlap with recognition sites.
        """
        indices = []  # No IUPAC codes.
        fwd_sites = [5, 15]  # Forward recognition sites.
        rev_sites = [25]  # Reverse recognition site.
        enzyme_oh_length = 6  # Enzyme overhang length.
        assert not check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # No overlap expected.

    def test_IUPAC_code_outside_sites(self):
        """
        Test case where IUPAC codes are outside the recognition sites.
        - The function should return False, as the IUPAC codes are not within the defined sites.
        """
        indices = [(3, 'R'), (20, 'Y')]  # IUPAC codes outside the recognition sites.
        fwd_sites = [5, 10]  # Forward recognition sites.
        rev_sites = [25]  # Reverse recognition site.
        enzyme_oh_length = 6  # Enzyme overhang length.
        assert not check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # No overlap expected.

    def test_IUPAC_code_in_fwd_site(self):
        """
        Test case where an IUPAC code overlaps a forward recognition site.
        - The function should return True, as the IUPAC code is inside the forward recognition site.
        """
        indices = [(6, 'R')]  # IUPAC code within the forward recognition site.
        fwd_sites = [5, 15]  # Forward recognition sites.
        rev_sites = [25]  # Reverse recognition site.
        enzyme_oh_length = 6  # Enzyme overhang length.
        assert check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # Overlap expected.

    def test_IUPAC_code_in_rev_site(self):
        """
        Test case where an IUPAC code overlaps a reverse recognition site.
        - The function should return True, as the IUPAC code is inside the reverse recognition site.
        """
        indices = [(27, 'Y')]  # IUPAC code within the reverse recognition site.
        fwd_sites = [5, 15]  # Forward recognition sites.
        rev_sites = [25]  # Reverse recognition site.
        enzyme_oh_length = 6  # Enzyme overhang length.
        assert check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # Overlap expected.

    def test_multiple_recognition_sites(self):
        """
        Test case with multiple recognition sites.
        - The function should return True if any IUPAC code overlaps with any recognition site.
        """
        indices = [(16, 'R'), (30, 'Y')]  # IUPAC codes within the recognition sites.
        fwd_sites = [5, 15]  # Forward recognition sites.
        rev_sites = [25, 35]  # Reverse recognition sites.
        enzyme_oh_length = 6  # Enzyme overhang length.
        assert check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # Overlap expected.

    def test_edge_case_overlap_end_site(self):
        """
        Test case where an IUPAC code is at the end of a recognition site.
        - The function should return True, as the IUPAC code overlaps the end of the recognition site.
        """
        indices = [(9, 'R')]  # IUPAC code at the end of a recognition site.
        fwd_sites = [5]  # Forward recognition site.
        rev_sites = [25]  # Reverse recognition site.
        enzyme_oh_length = 5  # Enzyme overhang length.
        assert check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # Overlap at the boundary.

    def test_edge_case_no_overlap_end_site(self):
        """
        Test case where an IUPAC code is at the end of a recognition site but does not overlap.
        - The function should return False, as the IUPAC code does not overlap the recognition site.
        """
        indices = [(10, 'R')]  # IUPAC code at the end, but no overlap.
        fwd_sites = [5]  # Forward recognition site.
        rev_sites = [25]  # Reverse recognition site.
        enzyme_oh_length = 5  # Enzyme overhang length.
        assert not check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)  # No overlap expected.


class TestCheckRecognitionSitesInExpandedSequences:
    @pytest.fixture
    def enzyme(self):
        """
        Fixture to create a test enzyme object.

        This fixture creates an enzyme object, which can be used in multiple tests.
        The enzyme is instantiated with a name, forward and reverse recognition sites,
        spacer length, and overhang length.

        Returns:
            Enzyme: A test enzyme object.
        """
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 0, 6)

    def test_no_IUPAC_code_overlap(self, enzyme):
        """
        Test when no IUPAC codes overlap with recognition sites.

        In this test, an enzyme is used to check if any recognition sites in the expanded
        sequences are affected by the IUPAC codes. Since no IUPAC codes overlap with the
        recognition sites, the expected result should be a list of False values.

        Args:
            enzyme (Enzyme): The enzyme object used for testing.

        Asserts:
            result (list): A list of booleans indicating whether recognition sites were affected.
        """
        enzyme = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        indices = [(5, 'N')]  # IUPAC code at index 5
        expanded_sequences = ["ATGCATGCAT", "GCTAGCTAAG"]

        result = check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme)

        # Assertions to check that no recognition sites are affected
        assert result == [False, False]

    def test_IUPAC_code_overlap(self, enzyme):
        """
        Test when an IUPAC code overlaps with a recognition site.

        In this test, the expanded sequences contain an IUPAC code that overlaps with
        the enzyme's recognition site. The expected result is that the recognition
        site will be identified as affected in both sequences.

        Args:
            enzyme (Enzyme): The enzyme object used for testing.

        Asserts:
            result (list): A list of booleans indicating whether recognition sites were affected.
        """
        enzyme = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        indices = [(5, 'N')]  # IUPAC code at index 5
        expanded_sequences = ["ATGCGAATTC", "GCTACTTAAG"]

        result = check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme)

        # Assertions to check that recognition sites are affected
        assert result == [True, True]

    def test_multiple_IUPAC_codes(self, enzyme):
        """
        Test with multiple IUPAC codes.

        This test checks whether the function can handle multiple IUPAC codes in
        a sequence. The enzyme's recognition site is checked against the positions
        of multiple IUPAC codes in the expanded sequences.

        Args:
            enzyme (Enzyme): The enzyme object used for testing.

        Asserts:
            result (list): A list of booleans indicating whether recognition sites were affected.
        """
        enzyme = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        indices = [(5, 'N'), (8, 'R')]  # Multiple IUPAC codes at indices 5 and 8
        expanded_sequences = ["ATGCGAATTC", "ATGCCGTAAG"]

        result = check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme)

        # Assertions to check that only the first sequence has the recognition site affected
        assert result == [True, False]


class TestAppendValidSequences:
    @pytest.fixture
    def setup_data(self):
        """
        Fixture to provide common test data.

        This fixture prepares and returns common data to be used across multiple tests.
        It includes a sequence name, an original sequence, and a list of expanded sequences.

        Returns:
            tuple: A tuple containing the name, sequence, and expanded sequences.
        """
        name = "TestSequence"
        seq = "ATGCTN"
        expanded_sequences = ["ATGCTA", "ATGCTG", "ATGCTC"]
        return name, seq, expanded_sequences

    def test_with_recognition_sites(self, setup_data):
        """
        Test appending expanded sequences when recognition sites exist.

        This test ensures that when recognition sites are present (True values in `rec_sites`),
        the expanded sequences are appended to `new_rows` with the correct information.

        Args:
            setup_data (tuple): The test data provided by the `setup_data` fixture.

        Asserts:
            new_rows (list): The list of appended rows should match the expected values.
        """
        name, seq, expanded_sequences = setup_data
        rec_sites = [True, True, True]
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions to check that the correct rows are added to `new_rows`
        assert len(new_rows) == 3
        for i, row in enumerate(new_rows):
            assert row['name'] == name
            assert row['DNA'] == seq
            assert row['valid_mixed_bases'] == expanded_sequences[i]

    def test_without_recognition_sites(self, setup_data):
        """
        Test appending the original sequence when no recognition sites exist.

        This test ensures that when no recognition sites are found (False values in `rec_sites`),
        only the original sequence is appended to `new_rows` with the expected information.

        Args:
            setup_data (tuple): The test data provided by the `setup_data` fixture.

        Asserts:
            new_rows (list): Only one row should be added, with the original sequence.
        """
        name, seq, expanded_sequences = setup_data
        rec_sites = [False, False, False]
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions to check that the only row added is the original sequence
        assert len(new_rows) == 1
        assert new_rows[0]['name'] == name
        assert new_rows[0]['DNA'] == seq
        assert new_rows[0]['valid_mixed_bases'] == seq

    def test_mixed_recognition_sites(self, setup_data):
        """
        Test mixed recognition site results.

        This test verifies that when recognition sites are mixed (True and False values in `rec_sites`),
        the correct sequences are appended based on the recognition site status.

        Args:
            setup_data (tuple): The test data provided by the `setup_data` fixture.

        Asserts:
            new_rows (list): The rows should be correctly populated based on mixed recognition sites.
        """
        name, seq, expanded_sequences = setup_data
        rec_sites = [True, False, True]
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions to ensure that the sequences are correctly handled based on recognition sites
        assert len(new_rows) == 3
        assert new_rows[0]['valid_mixed_bases'] == "ATGCTA"
        assert new_rows[1]['valid_mixed_bases'] == "ATGCTG"
        assert new_rows[2]['valid_mixed_bases'] == "ATGCTC"

    def test_empty_expanded_sequences(self):
        """
        Test with no expanded sequences provided.

        This test ensures that if no expanded sequences are given, the original sequence is still
        appended to `new_rows` as the only valid sequence.

        Asserts:
            new_rows (list): Only the original sequence should be added.
        """
        name = "TestSequence"
        seq = "ATGCTN"
        expanded_sequences = []
        rec_sites = []
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions to ensure that when there are no expanded sequences, only the original sequence is added
        assert len(new_rows) == 1
        assert new_rows[0]['name'] == name
        assert new_rows[0]['DNA'] == seq
        assert new_rows[0]['valid_mixed_bases'] == seq
