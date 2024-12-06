import pytest
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
from split_sites_utils import (
    find_starts_of_consecutive_indices,
    adjust_oligo_lengths,
    calculate_oligo_lengths,
    generate_cassettes,
    is_valid_overhang,
    hamming_distance,
    find_constant_indices
)# Replace `your_module` with the actual module name

class TestFindStartsOfConsecutiveIndices:
    """Test suite for the `find_starts_of_consecutive_indices` function."""

    def test_empty_conserved_indices(self):
        """Test with an empty conserved indices list."""
        conserved_indices = []
        consecutive_count = 4
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == []

    def test_no_consecutive_indices(self):
        """Test with no consecutive indices."""
        conserved_indices = [1, 3, 5, 7]
        consecutive_count = 3
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == []

    def test_single_consecutive_set(self):
        """Test with one set of consecutive indices."""
        conserved_indices = [1, 2, 3, 4, 6, 7]
        consecutive_count = 4
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == [0]

    def test_multiple_consecutive_sets(self):
        """Test with multiple sets of consecutive indices."""
        conserved_indices = [1, 2, 3, 4, 5, 7, 8, 9, 10]
        consecutive_count = 4
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == [0, 1, 5]

    def test_exact_length_match(self):
        """Test when the list length matches the consecutive count."""
        conserved_indices = [1, 2, 3, 4]
        consecutive_count = 4
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == [0]

    def test_consecutive_count_greater_than_list_length(self):
        """Test when the consecutive count is greater than the list length."""
        conserved_indices = [1, 2, 3]
        consecutive_count = 4
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == []

    def test_mixed_consecutive_and_non_consecutive(self):
        """Test with a mix of consecutive and non-consecutive indices."""
        conserved_indices = [1, 2, 3, 5, 6, 7, 9, 10, 11]
        consecutive_count = 3
        assert find_starts_of_consecutive_indices(conserved_indices, consecutive_count) == [0, 3, 6]

class TestAdjustOligoLengths:
        """Test suite for the `adjust_oligo_lengths` function."""

        def test_basic_adjustment(self):
            """Test a standard case with valid inputs."""
            min_oligo_len = 50
            max_oligo_len = 100
            recognition_site_len = 4
            spacer_bps = 3
            OH_len = 2
            assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (
            32, 82)

        def test_min_length_less_than_adjustment(self):
            """Test when the min_oligo_len becomes less than 1 after adjustment."""
            min_oligo_len = 10
            max_oligo_len = 100
            recognition_site_len = 10
            spacer_bps = 5
            OH_len = 5
            assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (
            1, 60)

        def test_zero_spacer_and_OH(self):
            """Test when spacer_bps and OH_len are zero."""
            min_oligo_len = 50
            max_oligo_len = 100
            recognition_site_len = 6
            spacer_bps = 0
            OH_len = 0
            assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (
            38, 88)

        def test_zero_recognition_site_len(self):
            """Test when recognition_site_len is zero."""
            min_oligo_len = 50
            max_oligo_len = 100
            recognition_site_len = 0
            spacer_bps = 2
            OH_len = 3
            assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (
            40, 90)

        def test_min_equals_max(self):
            """Test when min_oligo_len equals max_oligo_len."""
            min_oligo_len = 100
            max_oligo_len = 100
            recognition_site_len = 5
            spacer_bps = 2
            OH_len = 3
            assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (
            80, 80)

        def test_min_and_max_below_adjustment(self):
            """Test when both min_oligo_len and max_oligo_len are less than the adjustment."""
            min_oligo_len = 20
            max_oligo_len = 30
            recognition_site_len = 10
            spacer_bps = 5
            OH_len = 5
            try:
                adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len)
            except ValueError as e:
                assert str(e) == (
                    "The max oligo length (30) is incompatible with the selected enzyme. "
                    "It must be at least 40 to accommodate the enzyme properties "
                    "(recognition site: 10, spacer: 5, overhang: 5)."
                )
            else:
                assert False, "Expected ValueError to be raised for incompatible max oligo length."

        def test_exactly_equal_to_adjustment(self):
            """Test when max_oligo_len is exactly equal to the required adjustment."""
            min_oligo_len = 40
            max_oligo_len = 40
            recognition_site_len = 10
            spacer_bps = 5
            OH_len = 5
            assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (
            1, 1)

        def test_max_oligo_less_than_adjustment(self):
            """Test when max_oligo_len is less than the adjustment, raising an error."""
            min_oligo_len = 50
            max_oligo_len = 20
            recognition_site_len = 6
            spacer_bps = 5
            OH_len = 3
            try:
                adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len)
            except ValueError as e:
                assert str(e) == (
                    "The max oligo length (20) is incompatible with the selected enzyme. "
                    "It must be at least 28 to accommodate the enzyme properties "
                    "(recognition site: 6, spacer: 5, overhang: 3)."
                )
            else:
                assert False, "Expected ValueError to be raised for incompatible max oligo length."

class TestCalculateOligoLengths:
    """Test suite for the `calculate_oligo_lengths` function."""
    def test_no_split_indices(self):
        reference = "ATGCATGCATGC"
        split_indices = []
        OH_len = 5
        assert calculate_oligo_lengths(reference, split_indices, OH_len) == [12]

    def test_single_split_index(self):
        reference = "ATGCATGCATGC"
        split_indices = [4]
        OH_len = 5
        assert calculate_oligo_lengths(reference, split_indices, OH_len) == [9, 8]

    def test_multiple_split_indices(self):
        reference = "ATGCATGCATGC"
        split_indices = [4, 8]
        OH_len = 3
        assert calculate_oligo_lengths(reference, split_indices, OH_len) == [7, 7, 4]

    def test_invalid_split_index_negative(self):
        reference = "ATGCATGCATGC"
        split_indices = [-1, 4]
        OH_len = 2
        with pytest.raises(ValueError,
                           match="One or more split indices are out of the range of the reference sequence length."):
            calculate_oligo_lengths(reference, split_indices, OH_len)

    def test_invalid_split_index_out_of_range(self):
        reference = "ATGCATGCATGC"
        split_indices = [4, 13]
        OH_len = 2
        with pytest.raises(ValueError,
                           match="One or more split indices are out of the range of the reference sequence length."):
            calculate_oligo_lengths(reference, split_indices, OH_len)

    def test_split_at_start(self):
        reference = "ATGCATGCATGC"
        split_indices = [0]
        OH_len = 2
        assert calculate_oligo_lengths(reference, split_indices, OH_len) == [2, 12]

    def test_split_at_end(self):
        reference = "ATGCATGCATGC"
        split_indices = [12]
        OH_len = 4
        assert calculate_oligo_lengths(reference, split_indices, OH_len) == [16, 0]

class TestAdjustOligoLengths:
    """Test suite for the `adjust_oligo_lengths` function."""

    def test_basic_adjustment(self):
        """Test a standard case with valid inputs."""
        min_oligo_len = 50
        max_oligo_len = 100
        recognition_site_len = 4
        spacer_bps = 3
        OH_len = 2
        assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (32, 82)

    def test_min_length_less_than_adjustment(self):
        """Test when the min_oligo_len becomes less than 1 after adjustment."""
        min_oligo_len = 10
        max_oligo_len = 100
        recognition_site_len = 10
        spacer_bps = 5
        OH_len = 5
        assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (1, 60)

    def test_zero_spacer_and_OH(self):
        """Test when spacer_bps and OH_len are zero."""
        min_oligo_len = 50
        max_oligo_len = 100
        recognition_site_len = 6
        spacer_bps = 0
        OH_len = 0
        assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (38, 88)

    def test_zero_recognition_site_len(self):
        """Test when recognition_site_len is zero."""
        min_oligo_len = 50
        max_oligo_len = 100
        recognition_site_len = 0
        spacer_bps = 2
        OH_len = 3
        assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (40, 90)

    def test_min_equals_max(self):
        """Test when min_oligo_len equals max_oligo_len."""
        min_oligo_len = 100
        max_oligo_len = 100
        recognition_site_len = 5
        spacer_bps = 2
        OH_len = 3
        assert adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len) == (80, 80)

    def test_min_and_max_below_adjustment(self):
        """Test when both min_oligo_len and max_oligo_len are less than the adjustment."""
        min_oligo_len = 20
        max_oligo_len = 30
        recognition_site_len = 10
        spacer_bps = 5
        OH_len = 5
        try:
            adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len)
        except ValueError as e:
            assert str(e) == (
                "The max oligo length (30) is incompatible with the selected enzyme. "
                "It must be at least 41 to accommodate the enzyme properties "
                "(recognition site: 10, spacer: 5, overhang: 5)."
            )
        else:
            assert False, "Expected ValueError to be raised for incompatible max oligo length."

    def test_exactly_equal_to_adjustment(self):
        """Test when max_oligo_len is exactly equal to the required adjustment."""
        min_oligo_len = 40
        max_oligo_len = 40
        recognition_site_len = 10
        spacer_bps = 5
        OH_len = 5
        try:
            adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len)
        except ValueError as e:
            assert str(e) == (
                "The max oligo length (40) is incompatible with the selected enzyme. "
                "It must be at least 41 to accommodate the enzyme properties "
                "(recognition site: 10, spacer: 5, overhang: 5)."
            )
        else:
            assert False, "Expected ValueError to be raised for incompatible max oligo length."

    def test_max_oligo_len_less_than_adjustment(self):
        """Test when max_oligo_len is less than or equal to the adjustment, raising an error."""
        min_oligo_len = 50
        max_oligo_len = 20
        recognition_site_len = 6
        spacer_bps = 5
        OH_len = 3
        try:
            adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len)
        except ValueError as e:
            assert str(e) == (
                "The max oligo length (20) is incompatible with the selected enzyme. "
                "It must be at least 29 to accommodate the enzyme properties "
                "(recognition site: 6, spacer: 5, overhang: 3)."
            )
        else:
            assert False, "Expected ValueError to be raised for incompatible max oligo length."


class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = OH_length

    def __repr__(self):
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, OH_length={self.OH_length})"


# Test class for the generate_cassettes function
class TestGenerateCassettes:
    @pytest.fixture
    def enzyme(self):
        """Create a test enzyme."""
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 0, 3)

    def test_single_split_site(self, enzyme):
        DNA_df = pd.DataFrame({
            'rec_sites_removed': ['ATGCATGCATGC', 'GCATGCATGCAT']
        })
        split_sites = [0, 4, 12]

        # Call the generate_cassettes function with the enzyme's OH_length
        result = generate_cassettes(DNA_df, split_sites, enzyme)

        expected = pd.DataFrame({
            'rec_sites_removed': ['ATGCATGCATGC', 'GCATGCATGCAT'],
            'Cassette 1': ['ATGCATG', 'GCATGCA'],
            'Cassette 2': ['ATGCATGC', 'GCATGCAT'],
        })

        pd.testing.assert_frame_equal(result, expected)

    def test_multiple_split_sites(self, enzyme):
        DNA_df = pd.DataFrame({
            'rec_sites_removed': ['ATGCATGCATGC', 'GCATGCATGCAT'],
        })
        split_sites = [0, 2, 5, 12]

        # Call the generate_cassettes function with the enzyme's OH_length
        result = generate_cassettes(DNA_df, split_sites, enzyme)

        expected = pd.DataFrame({
            'rec_sites_removed': ['ATGCATGCATGC', 'GCATGCATGCAT'],
            'Cassette 1': ['ATGCA', 'GCATG'],
            'Cassette 2': ['GCATGC', 'ATGCAT'],
            'Cassette 3': ['TGCATGC', 'CATGCAT'],
        })

        pd.testing.assert_frame_equal(result, expected)

    def test_no_split_sites(self, enzyme):
        DNA_df = pd.DataFrame({
            'rec_sites_removed': ['ATGCATGC', 'GCATGCAT']
        })
        split_sites = [0,8]

        # Call the generate_cassettes function with the enzyme's OH_length
        result = generate_cassettes(DNA_df, split_sites, enzyme)

        expected = pd.DataFrame({
            'rec_sites_removed': ['ATGCATGC', 'GCATGCAT'],
            'Cassette 1': ['ATGCATGC', 'GCATGCAT']
        })

        pd.testing.assert_frame_equal(result, expected)


    def test_empty_dataframe(self, enzyme):
        DNA_df = pd.DataFrame({
            'rec_sites_removed': []
        })
        split_sites = [0, 4]

        # Call the generate_cassettes function with the enzyme's OH_length
        result = generate_cassettes(DNA_df, split_sites, enzyme)

        expected = pd.DataFrame({
            'rec_sites_removed': [],
            'Cassette 1': []
        })

        pd.testing.assert_frame_equal(result, expected)


class TestIsValidOverhang:

    def test_palindrome(self):
        """Test that palindromes return False"""
        overhang = "AGTGA"
        assert is_valid_overhang(overhang) is False

    def test_repeated_nucleotides(self):
        """Test that overhangs with three consecutive repeated nucleotides return False"""
        overhang = "AAAGT"
        assert is_valid_overhang(overhang) is False

    def test_extreme_gc_content_0(self):
        """Test that overhangs with 0% GC content return False"""
        overhang = "AAA"
        assert is_valid_overhang(overhang) is False

    def test_extreme_gc_content_100(self):
        """Test that overhangs with 100% GC content return False"""
        overhang = "GGG"
        assert is_valid_overhang(overhang) is False

    def test_valid_overhang(self):
        """Test that valid overhangs return True"""
        overhang = "AGCTG"
        assert is_valid_overhang(overhang) is True

    def test_edge_case_empty_overhang(self):
        """Test that an empty overhang returns False"""
        overhang = ""
        assert is_valid_overhang(overhang) is True

    def test_edge_case_single_nucleotide(self):
        """Test that a single nucleotide overhang is valid"""
        overhang = "A"
        assert is_valid_overhang(overhang) is True

    def test_mixed_gc_content(self):
        """Test that overhangs with mixed GC content return True"""
        overhang = "AGTCG"
        assert is_valid_overhang(overhang) is True


class TestHammingDistance:

    def test_equal_sequences(self):
        """Test if Hamming distance is 0 for identical sequences."""
        seq1 = "AGCTAG"
        seq2 = "AGCTAG"
        assert hamming_distance(seq1, seq2) == 0

    def test_completely_different_sequences(self):
        """Test if Hamming distance equals the length of the sequences when completely different."""
        seq1 = "AGCTAG"
        seq2 = "TCGATC"
        assert hamming_distance(seq1, seq2) == 6

    def test_one_difference(self):
        """Test Hamming distance for sequences that differ by one position."""
        seq1 = "AGCTAG"
        seq2 = "AGCTTG"
        assert hamming_distance(seq1, seq2) == 1

    def test_empty_sequences(self):
        """Test Hamming distance for empty sequences."""
        seq1 = ""
        seq2 = ""
        assert hamming_distance(seq1, seq2) == 0


class TestFindConstantIndices:

    def test_identical_sequences(self):
        """Test when all sequences are identical to the reference."""
        reference = "AGCTAG"
        sequences = ["AGCTAG", "AGCTAG", "AGCTAG"]
        expected = [0, 1, 2, 3, 4, 5]
        assert find_constant_indices(reference, sequences) == expected

    def test_partial_match(self):
        """Test when sequences have some matching positions with the reference."""
        reference = "AGCTAG"
        sequences = ["AGCTGG", "TGCTAG", "AGCGAG"]
        expected = [1, 2, 5]
        assert find_constant_indices(reference, sequences) == expected

    def test_no_match(self):
        """Test when no positions match between the reference and the sequences."""
        reference = "AGCTAG"
        sequences = ["TTTAAA", "GGGGGG"]
        expected = []
        assert find_constant_indices(reference, sequences) == expected

    def test_empty_sequences(self):
        """Test when sequences list is empty."""
        reference = "AGCTAG"
        sequences = []
        expected = [0, 1, 2, 3, 4, 5]
        assert find_constant_indices(reference, sequences) == expected

    def test_single_sequence(self):
        """Test when there is only one sequence in the list."""
        reference = "AGCTAG"
        sequences = ["AGCTGG"]
        expected = [0, 1, 2, 3, 5]
        assert find_constant_indices(reference, sequences) == expected

    def test_edge_case_empty_reference(self):
        """Test when the reference sequence is empty."""
        reference = ""
        sequences = ["", ""]
        expected = []
        assert find_constant_indices(reference, sequences) == expected

    def test_edge_case_one_char_reference(self):
        """Test when the reference sequence is one character long."""
        reference = "A"
        sequences = ["A", "A", "A"]
        expected = [0]
        assert find_constant_indices(reference, sequences) == expected
