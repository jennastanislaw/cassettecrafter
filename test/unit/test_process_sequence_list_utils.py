"""
Test script for validating sequence processing utilities used for DNA sequence
generation, modification, and validation.

Dependencies:
- pytest: Testing framework for unit tests.
- Bio.Seq (Biopython): Provides tools for manipulating DNA sequences, including
  the Seq class for representing DNA sequences.
- pandas: For handling DNA sequence data stored in DataFrame format.
- os, sys: For managing file paths and module imports.

Classes:
- TestGenerateRandomDNA: Contains tests for the generate_random_dna function.
- TestCheckRandomOligoForSites: Contains tests for the check_random_oligo_for_sites
  function.
- TestAddEndRestrictionSites: Contains tests for adding restriction enzyme recognition
  sites to DNA sequences in a DataFrame.
- TestGenerateSpacersRestSites: Contains tests for generating DNA sequences with
  restriction enzyme recognition sites and spacers.
- TestGenerateSpacer: Contains tests for generating spacer sequences.
- TestGenerateSitelessSequenceIntegration: Contains tests for ensuring recognition
  sites are removed from DNA sequences.

Note:
- The Enzyme class is used to model restriction enzymes, including their forward
  and reverse recognition sites, spacer lengths, and other relevant properties.
- Test functions are designed to check various edge cases, such as empty DNA sequences,
  sequences without recognition sites, and sequences that contain multiple recognition sites.
"""

import pytest
from Bio.Seq import Seq
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from process_sequence_list_utils import (
    generate_random_dna,
    check_random_oligo_for_sites,
    add_end_restriction_sites,
    generate_spacers_rest_sites,
    generate_spacer,
    generate_siteless_sequence,
    ensure_minimum_length
)


class TestGenerateRandomDNA:
    def test_generate_random_dna_valid_length(self):
        """Test generation of a valid-length DNA sequence."""
        length = 10
        dna_sequence = generate_random_dna(length)
        # Ensure sequence is a Biopython Seq object
        assert isinstance(dna_sequence, Seq), "The result should be a Biopython Seq object."
        # Verify sequence length matches the requested length
        assert len(dna_sequence) == length, f"Expected length {length}, got {len(dna_sequence)}."
        # Ensure sequence contains only valid bases (A, T, C, G)
        assert all(base in "ATCG" for base in dna_sequence), "Sequence should only contain A, T, C, G."

    def test_generate_random_dna_zero_length(self):
        """Test generation of an empty DNA sequence for zero length."""
        dna_sequence = generate_random_dna(0)
        # Ensure result is a Biopython Seq object
        assert isinstance(dna_sequence, Seq), "The result should be a Biopython Seq object."
        # Verify sequence is empty for length zero
        assert len(dna_sequence) == 0, "The sequence should be empty for length zero."

    def test_generate_random_dna_negative_length(self):
        """Test that negative length raises an error."""
        with pytest.raises(ValueError, match="Length of DNA sequence must be a positive integer."):
            generate_random_dna(-5)

    def test_generate_random_dna_randomness(self):
        """Test that generated sequences are random."""
        length = 10
        sequence1 = generate_random_dna(length)
        sequence2 = generate_random_dna(length)
        # Ensure the generated sequences are not identical
        assert sequence1 != sequence2, "Generated sequences should not be identical."


class Enzyme:
    """
    Class representing a restriction enzyme with its properties.

    Attributes:
        name (str): Name of the enzyme.
        fwd_recognition_site (str): Forward recognition site sequence.
        rev_recognition_site (str): Reverse recognition site sequence.
        spacer_length (int): Length of the spacer sequence.
        OH_length (int): Length of the overhang sequence.
    """

    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = OH_length

    def __repr__(self):
        """Return a string representation of the enzyme."""
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, OH_length={self.OH_length})"


class TestCheckRandomOligoForSites:
    @pytest.fixture
    def enzyme(self):
        """Create and return a test enzyme object."""
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 0, 4)

    def test_check_random_oligo_with_no_sites(self, enzyme):
        """Test that DNA sequence with no recognition sites returns 0 matches."""
        dna_sequence = "TTGACGTTGACGTTGACG"  # No recognition site
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 0, f"Expected 0 matches, got {result}."

    def test_check_random_oligo_with_one_site(self, enzyme):
        """Test that DNA sequence with one recognition site returns 1 match."""
        dna_sequence = "TTGACGAATTCGACGTTGA"  # One forward recognition site
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 1, f"Expected 1 match, got {result}."

    def test_check_random_oligo_with_multiple_sites(self, enzyme):
        """Test that DNA sequence with multiple recognition sites returns correct match count."""
        dna_sequence = "GAATTCCTTAAGGAATTC"  # Two forward and one reverse recognition sites
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 3, f"Expected 3 matches, got {result}."

    def test_empty_dna_sequence(self, enzyme):
        """Test that an empty DNA sequence returns 0 matches."""
        dna_sequence = ""
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 0, "Expected 0 matches for an empty DNA sequence."

    def test_check_random_oligo_with_only_reverse_sites(self, enzyme):
        """Test that DNA sequence with only reverse recognition sites returns correct match count."""
        dna_sequence = "CTTAAGCTTAAG"  # Two reverse recognition sites
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 2, f"Expected 2 matches, got {result}."


class TestAddEndRestrictionSites:
    """
        Test suite for the function `add_end_restriction_sites` which adds restriction enzyme recognition sites
        and spacer sequences to DNA sequences in a DataFrame.

        This test class checks the following functionalities:
        - Correct addition of forward and reverse recognition sites.
        - Correct spacer length for the enzyme.
        - Ensuring all sequences are processed and updated correctly in the DataFrame.

        Fixtures:
            enzyme: Creates a test enzyme object with recognition sites and spacer length.
            dna_dataframe: Creates a sample DataFrame with DNA sequences to test the function.

        Tests:
            - test_add_restriction_sites_integration: Validates that the forward and reverse recognition sites
              are correctly added to DNA sequences.
            - test_length: Ensures that the spacer length is correctly applied when adding restriction sites.
            - test_all_sequences_processed: Verifies that all DNA sequences are processed and updated correctly.
        """

    @pytest.fixture
    def enzyme(self):
        """Create and return a test enzyme object."""
        # Enzyme with specific recognition sites and spacer length
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 4, 6)

    @pytest.fixture
    def dna_dataframe(self):
        """Create a sample DataFrame for testing with DNA sequences."""
        # Sample DataFrame with DNA sequences for testing
        data = {
            "name": ["Sequence1", "Sequence2", "Sequence3"],
            "DNA": ["ATGCGTACG", "TGCATGCA", "CGTACGTAGC"],
            "DNA_unchanged": ["ATGCGTACG", "TGCATGCA", "CGTACGTAGC"]
        }
        return pd.DataFrame(data)

    def test_add_restriction_sites_integration(self, dna_dataframe, enzyme):
        """Test adding restriction sites to DNA sequences, ensuring recognition sites are added."""
        # Apply the function to the DataFrame
        result_df = add_end_restriction_sites(dna_dataframe, "DNA", enzyme)

        # Check if the forward and reverse recognition sites are present in each DNA sequence
        for index, row in result_df.iterrows():
            dna_sequence = row["DNA"]
            assert enzyme.fwd_recognition_site in dna_sequence, f"Forward recognition site missing in {dna_sequence}"
            assert enzyme.rev_recognition_site in dna_sequence, f"Reverse recognition site missing in {dna_sequence}"

    def test_length(self, dna_dataframe, enzyme):
        """Test that the spacer length is correctly applied when adding restriction sites."""
        # Apply the function to the DataFrame
        result_df = add_end_restriction_sites(dna_dataframe, "DNA", enzyme)

        # Check if the DNA sequence length matches the expected length, considering spacer and recognition site lengths
        for index, row in result_df.iterrows():
            dna_sequence = row["DNA"]
            expected_length = enzyme.spacer_length * 2 + len(enzyme.rev_recognition_site) + len(
                enzyme.fwd_recognition_site) + len(row["DNA_unchanged"])
            # Assert that the DNA sequence length matches the expected length
            assert len(dna_sequence) == expected_length, f"Length mismatch in {dna_sequence}"

    def test_all_sequences_processed(self, dna_dataframe, enzyme):
        """Test that all DNA sequences are processed and updated correctly in the DataFrame."""
        # Apply the function to the DataFrame
        result_df = add_end_restriction_sites(dna_dataframe, "DNA", enzyme)

        # Ensure the number of sequences remains the same in the result DataFrame
        assert len(result_df) == len(dna_dataframe), "Mismatch in the number of sequences"

        # Ensure each sequence has been updated correctly (by checking against the unchanged sequences)
        for index, row in result_df.iterrows():
            assert row["DNA"] != dna_dataframe["DNA_unchanged"].iloc[
                index], f"Sequence at index {index} was not updated correctly"


class TestGenerateSpacersRestSites:
    """
        Test suite for the function `generate_spacers_rest_sites` which generates DNA sequences with added
        restriction enzyme recognition sites and spacer sequences for a given enzyme.

        This test class checks the following functionalities:
        - Correct structure of the generated sequence with recognition sites and spacer sequences.
        - Correct spacer lengths for both forward and reverse spacers.
        - Handling of empty DNA sequences by ensuring the correct placement of recognition sites and spacers.

        Fixtures:
            enzyme: Creates a test enzyme object with recognition sites and spacer length.

        Tests:
            - test_valid_dna_sequence: Verifies that the function correctly adds recognition sites and spacer sequences
              to a valid DNA sequence.
            - test_empty_dna_sequence: Verifies that the function handles an empty DNA sequence properly.
        """

    @pytest.fixture
    def enzyme(self):
        """Create and return a test enzyme object with spacer lengths and recognition sites."""
        # Create an enzyme object with specific recognition sites, spacer length, and OH length
        return Enzyme(
            name="TestEnzyme",
            fwd_recognition_site="GAATTC",
            rev_recognition_site="CTTAAG",
            spacer_length=4,
            OH_length=6
        )

    def test_valid_dna_sequence(self, enzyme):
        """Test the function with a valid DNA sequence."""
        # Test DNA sequence
        original_dna = "ATGCGTACG"

        # Generate the sequence with spacer and restriction sites
        result = generate_spacers_rest_sites(original_dna, enzyme)

        # Verify the structure of the result: correct placement of recognition sites and spacers
        assert result.startswith(enzyme.fwd_recognition_site), "Forward recognition site missing"
        assert result.endswith(enzyme.rev_recognition_site), "Reverse recognition site missing"

        # Check that the forward and reverse spacer lengths are correct
        forward_spacer = result[
                         len(enzyme.fwd_recognition_site):len(enzyme.fwd_recognition_site) + enzyme.spacer_length]
        reverse_spacer = result[
                         -(len(enzyme.rev_recognition_site) + enzyme.spacer_length):-len(enzyme.rev_recognition_site)]
        assert len(forward_spacer) == enzyme.spacer_length, "Forward spacer length is incorrect"
        assert len(reverse_spacer) == enzyme.spacer_length, "Reverse spacer length is incorrect"

        # Check that the middle DNA sequence is placed correctly
        middle_dna = result[len(enzyme.fwd_recognition_site) + enzyme.spacer_length:-(
                len(enzyme.rev_recognition_site) + enzyme.spacer_length)]
        assert middle_dna == original_dna, "Original DNA sequence is not placed correctly"

    def test_empty_dna_sequence(self, enzyme):
        """Test the function with an empty DNA sequence."""
        original_dna = ""

        # Generate the sequence with spacer and restriction sites (should be empty DNA)
        result = generate_spacers_rest_sites(original_dna, enzyme)

        # Verify the structure of the result: correct placement of recognition sites and spacers
        assert result.startswith(enzyme.fwd_recognition_site), "Forward recognition site missing"
        assert result.endswith(enzyme.rev_recognition_site), "Reverse recognition site missing"

        # Check that spacer lengths are correct (even with empty DNA)
        forward_spacer = result[
                         len(enzyme.fwd_recognition_site):len(enzyme.fwd_recognition_site) + enzyme.spacer_length]
        reverse_spacer = result[
                         -(len(enzyme.rev_recognition_site) + enzyme.spacer_length):-len(enzyme.rev_recognition_site)]
        assert len(forward_spacer) == enzyme.spacer_length, "Forward spacer length is incorrect"
        assert len(reverse_spacer) == enzyme.spacer_length, "Reverse spacer length is incorrect"

        # Ensure no additional DNA is included in the result
        middle_dna = result[len(enzyme.fwd_recognition_site) + enzyme.spacer_length:-(
                len(enzyme.rev_recognition_site) + enzyme.spacer_length)]
        assert middle_dna == "", "Unexpected DNA sequence found"


class TestGenerateSpacer:
    """
    Test suite for the function `generate_spacer` which generates spacer sequences in DNA based on given enzyme
    recognition sites and directions.

    This test class checks the following functionalities:
    - Correct generation of forward and reverse spacers.
    - Validation of the complete DNA sequence with the spacer and recognition sites.
    - Error handling for invalid direction inputs.
    - End-to-end testing using random DNA sequences.

    Fixtures:
        enzyme: Creates a test enzyme object with recognition sites and spacer length.

    Tests:
        - test_generate_spacer_with_valid_inputs: Verifies that forward spacers are correctly generated and validated.
        - test_generate_reverse_spacer_with_valid_inputs: Verifies that reverse spacers are correctly generated and validated.
        - test_generate_spacer_handles_invalid_direction: Ensures that an error is raised for invalid direction input.
        - test_end_to_end_with_random_dna: Tests the full workflow using randomly generated DNA sequences.
    """

    @pytest.fixture
    def enzyme(self):
        """Fixture to create a test enzyme."""
        return Enzyme(
            name="TestEnzyme",
            fwd_recognition_site="GAATTC",
            rev_recognition_site="CTTAAG",
            spacer_length=4,
            OH_length=6  # Corrected attribute name
        )

    def test_generate_spacer_with_valid_inputs(self, enzyme):
        """Test that the spacer is correctly generated and validated."""
        dna = "ATGCGTACG"
        rec_site = enzyme.fwd_recognition_site

        # Generate a spacer in the forward direction
        spacer = generate_spacer(dna, rec_site, enzyme, direction="forward")

        # Verify the complete DNA sequence
        dna_with_spacer = rec_site + spacer + dna
        assert len(spacer) == enzyme.spacer_length, "Spacer length mismatch"
        assert check_random_oligo_for_sites(dna_with_spacer,
                                            enzyme) == 1, "Recognition site check failed for forward spacer"

    def test_generate_reverse_spacer_with_valid_inputs(self, enzyme):
        """Test that reverse spacers are correctly generated and validated."""
        dna = "ATGCGTACG"
        rec_site = enzyme.rev_recognition_site

        # Generate a spacer in the reverse direction
        spacer = generate_spacer(dna, rec_site, enzyme, direction="reverse")

        # Verify the complete DNA sequence
        dna_with_spacer = dna + spacer + rec_site
        assert len(spacer) == enzyme.spacer_length, "Spacer length mismatch"
        assert check_random_oligo_for_sites(dna_with_spacer,
                                            enzyme) == 1, "Recognition site check failed for reverse spacer"

    def test_generate_spacer_handles_invalid_direction(self, enzyme):
        """Test that invalid direction raises an error."""
        dna = "ATGCGTACG"
        rec_site = enzyme.fwd_recognition_site

        with pytest.raises(ValueError, match="Invalid direction"):
            generate_spacer(dna, rec_site, enzyme, direction="sideways")

    def test_end_to_end_with_random_dna(self, enzyme):
        """Test full workflow using a randomly generated DNA sequence."""
        random_dna = generate_random_dna(20)  # Generate random DNA
        forward_spacer = generate_spacer(random_dna, enzyme.fwd_recognition_site, enzyme, direction="forward")
        reverse_spacer = generate_spacer(random_dna, enzyme.rev_recognition_site, enzyme, direction="reverse")

        # Combine into a complete sequence
        full_sequence = enzyme.fwd_recognition_site + forward_spacer + random_dna + reverse_spacer + enzyme.rev_recognition_site
        assert check_random_oligo_for_sites(full_sequence, enzyme) == 2, "End-to-end spacer generation failed"


class TestGenerateSitelessSequenceIntegration:
    """
    Test suite for the function `generate_siteless_sequence` which generates random DNA sequences that do not
    contain any recognition sites for the given enzyme.

    This test class checks the following functionalities:
    - Correct generation of siteless DNA sequences without any recognition sites.
    - Ensuring the sequence length matches the desired length.
    - Ensuring randomness of generated sequences.
    - Correct handling of edge cases such as minimum sequence length.

    Fixtures:
        enzyme: Creates a test enzyme object with recognition sites and spacer length.

    Tests:
        - test_generate_siteless_sequence_valid_length: Verifies that the generated sequence has the correct length.
        - test_generate_siteless_sequence_randomness: Ensures that multiple calls generate unique sequences.
        - test_generate_siteless_sequence_edge_case_min_length: Verifies that the function handles the shortest possible sequence length.
        - test_generate_siteless_sequence_handles_recognition_sites: Ensures that generated sequences do not contain enzyme recognition sites.
    """

    @pytest.fixture
    def enzyme(self):
        """Create a test enzyme."""
        return Enzyme(
            name="TestEnzyme",
            fwd_recognition_site="GAATTC",
            rev_recognition_site="CTTAAG",
            spacer_length=4,
            OH_length=6
        )

    def test_generate_siteless_sequence_valid_length(self, enzyme):
        """Test that the generated sequence has the correct length."""
        n = 20  # Desired length of the DNA sequence
        result = generate_siteless_sequence(n, enzyme)

        assert len(result) == n, f"Generated sequence length {len(result)} does not match expected length {n}"
        assert check_random_oligo_for_sites(result, enzyme) == 0, "Generated sequence contains enzyme recognition sites"

    def test_generate_siteless_sequence_randomness(self, enzyme):
        """Test that multiple calls produce unique sequences."""
        n = 20
        sequences = {generate_siteless_sequence(n, enzyme) for _ in range(10)}

        # Ensure all sequences are unique
        assert len(sequences) == 10, "Generated sequences are not unique"

    def test_generate_siteless_sequence_edge_case_min_length(self, enzyme):
        """Test that the function works with a very short sequence."""
        n = len(enzyme.fwd_recognition_site)  # Minimum sequence length equal to the recognition site
        result = generate_siteless_sequence(n, enzyme)

        assert len(result) == n, f"Generated sequence length {len(result)} does not match expected length {n}"
        assert check_random_oligo_for_sites(result, enzyme) == 0, "Generated sequence contains enzyme recognition sites"

    def test_generate_siteless_sequence_handles_recognition_sites(self, enzyme):
        """Test that the function filters out sequences with recognition sites."""
        n = 30
        for _ in range(10):
            result = generate_siteless_sequence(n, enzyme)

            # Ensure that the forward and reverse recognition sites are not in the generated sequence
            assert enzyme.fwd_recognition_site not in result, "Generated sequence contains the forward recognition site"
            assert enzyme.rev_recognition_site not in result, "Generated sequence contains the reverse recognition site"


class TestEnsureMinimumLength:
    """
        Tests for ensuring that DNA sequences in a DataFrame are modified to meet a specified minimum length.

        This class contains tests to verify the behavior of the `ensure_minimum_length` function when applied to
        sequences of varying lengths. It includes checks for the following cases:
        - Handling non-string sequences (e.g., list, tuple, integer) and converting them to strings.
        - Handling edge cases where sequences are exactly the minimum length.
        - Handling mixed sequences, where some are shorter and need modification, and others are long enough to remain unchanged.

        Attributes:
            enzyme (Enzyme): A test enzyme with forward and reverse recognition sites.
        """

    @pytest.fixture
    def enzyme(self):
        """Create a test enzyme."""
        return Enzyme(
            name="TestEnzyme",
            fwd_recognition_site="GAATTC",
            rev_recognition_site="CTTAAG",
            spacer_length=4,
            OH_length=6
        )

    def test_non_string_sequences(self, enzyme):
        """Test with non-string DNA sequences (e.g., list, tuple, int)."""
        df = pd.DataFrame({
            'Modified DNA': [[1, 2, 3], (4, 5, 6), 12345]  # Non-string sequences (list, tuple, int)
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        # Assert that all sequences are converted to strings and modified to meet the minimum length
        assert all(
            isinstance(seq, str) for seq in result['Modified DNA']), "Not all sequences were converted to strings"
        assert all(len(seq) >= min_oligo_size for seq in
                   result['Modified DNA']), f"Some sequences are shorter than {min_oligo_size}"

    def test_edge_case_min_length(self, enzyme):
        """Test with a sequence that is exactly equal to the minimum length."""
        df = pd.DataFrame({
            'Modified DNA': ['ATCGGT']  # Sequence is exactly 6 bases long
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        # Assert that the sequence is expanded to meet the minimum length
        assert len(
            result['Modified DNA'][0]) == 18, f"Sequence length {len(result['Modified DNA'][0])} is not as expected"

    def test_mixed_sequences(self, enzyme):
        """Test with a mix of sequences that need modification and those that don't."""
        df = pd.DataFrame({
            'Modified DNA': ['ATCG', 'AATTGCGC', 'A']  # Mixed lengths
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        # Assert that sequences that were too short are modified
        assert len(result['Modified DNA'][0]) >= min_oligo_size, f"Sequence 0 is shorter than {min_oligo_size}"
        assert len(result['Modified DNA'][1]) == 20, f"Sequence 1 was unexpectedly modified"
        assert len(result['Modified DNA'][2]) >= min_oligo_size, f"Sequence 2 is shorter than {min_oligo_size}"


