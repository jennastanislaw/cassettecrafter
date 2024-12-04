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
import random

class TestGenerateRandomDNA:
    def test_generate_random_dna_valid_length(self):
        """Test generation of a DNA sequence of valid length."""
        length = 10
        dna_sequence = generate_random_dna(length)
        assert isinstance(dna_sequence, Seq), "The result should be a Biopython Seq object."
        assert len(dna_sequence) == length, f"Expected length {length}, got {len(dna_sequence)}."
        assert all(base in "ATCG" for base in dna_sequence), "Sequence should only contain A, T, C, G."

    def test_generate_random_dna_zero_length(self):
        """Test generation of a DNA sequence with length zero returns an empty Seq object."""
        dna_sequence = generate_random_dna(0)
        assert isinstance(dna_sequence, Seq), "The result should be a Biopython Seq object."
        assert len(dna_sequence) == 0, "The sequence should be empty for length zero."

    def test_generate_random_dna_negative_length(self):
        """Test generation of a DNA sequence with a negative length raises an error."""
        with pytest.raises(ValueError, match="Length of DNA sequence must be a positive integer."):
            generate_random_dna(-5)

    def test_generate_random_dna_randomness(self):
        """Test that the function generates random sequences."""
        length = 10
        sequence1 = generate_random_dna(length)
        sequence2 = generate_random_dna(length)
        assert sequence1 != sequence2, "Generated sequences should not be identical."


class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = OH_length

    def __repr__(self):
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, OH_length={self.OH_length})"



class TestCheckRandomOligoForSites:
    @pytest.fixture
    def enzyme(self):
        """Create a test enzyme."""
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 0, 6)

    def test_check_random_oligo_with_no_sites(self, enzyme):
        """Test when the DNA has no recognition sites."""
        dna_sequence = "TTGACGTTGACGTTGACG"  # No recognition site
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 0, f"Expected 0 matches, got {result}."

    def test_check_random_oligo_with_one_site(self, enzyme):
        """Test when the DNA has one recognition site."""
        dna_sequence = "TTGACGAATTCGACGTTGA"  # One forward recognition site
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 1, f"Expected 1 match, got {result}."

    def test_check_random_oligo_with_multiple_sites(self, enzyme):
        """Test when the DNA has multiple recognition sites."""
        dna_sequence = "GAATTCCTTAAGGAATTC"  # Two forward and one reverse recognition sites
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 3, f"Expected 3 matches, got {result}."

    def test_empty_dna_sequence(self, enzyme):
        """Test with an empty DNA sequence."""
        dna_sequence = ""
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 0, "Expected 0 matches for an empty DNA sequence."

    def test_check_random_oligo_with_only_reverse_sites(self, enzyme):
        """Test when the DNA has only reverse recognition sites."""
        dna_sequence = "CTTAAGCTTAAG"  # Two reverse recognition sites
        result = check_random_oligo_for_sites(dna_sequence, enzyme)
        assert result == 2, f"Expected 2 matches, got {result}."

class TestAddEndRestrictionSites:
    @pytest.fixture
    def enzyme(self):
        """Create a test enzyme."""
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 4, 6)

    @pytest.fixture
    def dna_dataframe(self):
        """Create a sample DataFrame for testing."""
        data = {
            "name": ["Sequence1", "Sequence2", "Sequence3"],
            "DNA": ["ATGCGTACG", "TGCATGCA", "CGTACGTAGC"],
            "DNA_unchanged": ["ATGCGTACG", "TGCATGCA", "CGTACGTAGC"]
        }
        return pd.DataFrame(data)

    def test_add_restriction_sites_integration(self, dna_dataframe, enzyme):
        """Integration test for adding restriction sites and spacer sequences."""
        # Apply the function to the DataFrame
        result_df = add_end_restriction_sites(dna_dataframe, "DNA", enzyme)

        # Check if the recognition sites are correctly added
        for index, row in result_df.iterrows():
            dna_sequence = row["DNA"]
            assert enzyme.fwd_recognition_site in dna_sequence, f"Forward recognition site missing in {dna_sequence}"
            assert enzyme.rev_recognition_site in dna_sequence, f"Reverse recognition site missing in {dna_sequence}"

    def test_length(self, dna_dataframe, enzyme):
        """Test that the spacer length is correctly respected when adding restriction sites."""
        # Apply the function to the DataFrame
        result_df = add_end_restriction_sites(dna_dataframe, "DNA", enzyme)

        # Check that the spacer length (defined by the enzyme) is correct
        for index, row in result_df.iterrows():
            dna_sequence = row["DNA"]

            expected_length = enzyme.spacer_length * 2 + len(enzyme.rev_recognition_site) + len(
                enzyme.fwd_recognition_site) + len(row["DNA_unchanged"])
            # Calculate the expected spacer length by subtracting the lengths of the recognition sites
            assert len(dna_sequence) == expected_length, f"Length mismatch in {dna_sequence}"


    def test_all_sequences_processed(self, dna_dataframe, enzyme):
        """Test that all sequences are processed and updated correctly."""
        # Apply the function to the DataFrame
        result_df = add_end_restriction_sites(dna_dataframe, "DNA", enzyme)

        # Check that the output DataFrame has the same number of rows as the input DataFrame
        assert len(result_df) == len(dna_dataframe), "Mismatch in the number of sequences"

        # Ensure that all DNA sequences are updated
        for index, row in result_df.iterrows():
            assert row["DNA"] != dna_dataframe["DNA_unchanged"].iloc[index], f"Sequence at index {index} was not updated correctly"


class TestGenerateSpacersRestSites:
    @pytest.fixture
    def enzyme(self):
        """Fixture to create a test enzyme."""
        return Enzyme(
            name="TestEnzyme",
            fwd_recognition_site="GAATTC",
            rev_recognition_site="CTTAAG",
            spacer_length=4,
            OH_length=6
        )

    def test_valid_dna_sequence(self, enzyme):
        """Test the function with a valid DNA sequence."""
        original_dna = "ATGCGTACG"
        result = generate_spacers_rest_sites(original_dna, enzyme)

        # Verify the structure of the result
        assert result.startswith(enzyme.fwd_recognition_site), "Forward recognition site missing"
        assert result.endswith(enzyme.rev_recognition_site), "Reverse recognition site missing"

        # Check spacer lengths
        forward_spacer = result[len(enzyme.fwd_recognition_site):len(enzyme.fwd_recognition_site) + enzyme.spacer_length]
        reverse_spacer = result[-(len(enzyme.rev_recognition_site) + enzyme.spacer_length):-len(enzyme.rev_recognition_site)]
        assert len(forward_spacer) == enzyme.spacer_length, "Forward spacer length is incorrect"
        assert len(reverse_spacer) == enzyme.spacer_length, "Reverse spacer length is incorrect"

        # Check the DNA placement
        middle_dna = result[len(enzyme.fwd_recognition_site) + enzyme.spacer_length:-(len(enzyme.rev_recognition_site) + enzyme.spacer_length)]
        assert middle_dna == original_dna, "Original DNA sequence is not placed correctly"

    def test_empty_dna_sequence(self, enzyme):
        """Test the function with an empty DNA sequence."""
        original_dna = ""
        result = generate_spacers_rest_sites(original_dna, enzyme)

        # Verify the structure of the result
        assert result.startswith(enzyme.fwd_recognition_site), "Forward recognition site missing"
        assert result.endswith(enzyme.rev_recognition_site), "Reverse recognition site missing"

        # Check spacer lengths
        forward_spacer = result[len(enzyme.fwd_recognition_site):len(enzyme.fwd_recognition_site) + enzyme.spacer_length]
        reverse_spacer = result[-(len(enzyme.rev_recognition_site) + enzyme.spacer_length):-len(enzyme.rev_recognition_site)]
        assert len(forward_spacer) == enzyme.spacer_length, "Forward spacer length is incorrect"
        assert len(reverse_spacer) == enzyme.spacer_length, "Reverse spacer length is incorrect"

        # Ensure no additional DNA is included
        middle_dna = result[len(enzyme.fwd_recognition_site) + enzyme.spacer_length:-(len(enzyme.rev_recognition_site) + enzyme.spacer_length)]
        assert middle_dna == "", "Unexpected DNA sequence found"


class TestGenerateSpacer:
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

    def test_some_sequences_shorter(self, enzyme):
        """Test where some sequences are shorter than the minimum length"""
        df = pd.DataFrame({
            'Modified DNA': ['ATCG', 'AATTGCGC', 'A']  # First and third are shorter than 6
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        # Assert that the sequences that were too short are modified to meet the minimum length
        assert len(result['Modified DNA'][0]) >= min_oligo_size
        assert len(result['Modified DNA'][1]) == 20  # No change for long enough sequence
        assert len(result['Modified DNA'][2]) >= min_oligo_size

    def test_non_string_sequences(self, enzyme):
        """Test with non-string DNA sequences"""
        df = pd.DataFrame({
            'Modified DNA': [[1, 2, 3], (4, 5, 6), 12345]  # Non-string sequences (list, tuple, int)
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        # Assert that all sequences are converted to strings and modified to meet the minimum length
        assert isinstance(result['Modified DNA'][0], str)
        assert isinstance(result['Modified DNA'][1], str)
        assert isinstance(result['Modified DNA'][2], str)
        assert len(result['Modified DNA'][0]) >= min_oligo_size
        assert len(result['Modified DNA'][1]) >= min_oligo_size
        assert len(result['Modified DNA'][2]) >= min_oligo_size

    def test_edge_case_min_length(self, enzyme):
        """Test with a sequence that is exactly equal to the minimum length"""
        df = pd.DataFrame({
            'Modified DNA': ['ATCGGT']  # Sequence is exactly 6 bases long
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        assert len(result['Modified DNA'][0])== 18

    def test_mixed_sequences(self, enzyme):
        """Test with a mix of sequences that need modification and those that don't"""
        df = pd.DataFrame({
            'Modified DNA': ['ATCG', 'AATTGCGC', 'A']  # Mixed lengths
        })
        min_oligo_size = 6

        result = ensure_minimum_length(df, 'Modified DNA', min_oligo_size, enzyme)

        # Assert that sequences that were too short are modified
        assert len(result['Modified DNA'][0]) >= min_oligo_size
        assert len(result['Modified DNA'][1]) == 20  # No change for long enough sequence
        assert len(result['Modified DNA'][2]) >= min_oligo_size
