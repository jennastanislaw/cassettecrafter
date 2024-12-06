import pytest
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from split_sites import (
    find_split_indices,
    find_valid_overhangs,
    find_valid_combinations,
    calculate_hamming_distances
)

# Sample enzyme class for testing
class Enzyme:
    def __init__(self, fwd_recognition_site, OH_length, spacer_length):
        self.fwd_recognition_site = fwd_recognition_site
        self.OH_length = OH_length
        self.spacer_length = spacer_length

@pytest.fixture
def test_data():
    reference = "ATGCGTACGTTAGCTAGCGTACGTTAGC"
    sequences = [
        "ATGCGTACGTTAGCTAGCGTACGTTAGC",
        "ATGCGTACGTTAGCTAGCGTACGTTTGC",
        "ATGCGTACGTTAGCTAGCGTACGTTTAC"
    ]
    enzyme = Enzyme(fwd_recognition_site="GAATTC", OH_length=4, spacer_length=2)
    min_oligo_len = 6
    max_oligo_len = 37
    return reference, sequences, enzyme, min_oligo_len, max_oligo_len


class TestFindSplitIndices:
    def test_valid_output(self, test_data):
        reference, sequences, enzyme, min_oligo_len, max_oligo_len = test_data

        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # Verify the result has the expected format
        assert isinstance(result, list), "Output should be a list"
        assert all(isinstance(i, int) for i in result), "All indices should be integers"
        assert result[0] == 0, "First index should be 0"
        assert result[-1] == len(reference), "Last index should be the length of the reference"

    def test_respects_constraints(self, test_data):
        reference, sequences, enzyme, min_oligo_len, max_oligo_len = test_data

        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)
        OH_len = enzyme.OH_length

        result = result[1:-1]

        # Initialize the list for oligo lengths
        oligo_lengths = []

        # Calculate oligo lengths based on split indices
        if not result:
            # No splits, return the full length of the reference
            oligo_lengths.append(len(reference))
        else:
            # Calculate oligo lengths normally
            oligo_lengths.append(result[0] + OH_len)  # From start to first index
            oligo_lengths.extend(
                result[i + 1] - result[i] + OH_len for i in range(len(result) - 2)
            )
            oligo_lengths.append(len(reference) - result[-1])  # From last index to end

        # Check that all oligo lengths are within the specified range
        assert all(min_oligo_len <= length <= max_oligo_len for length in oligo_lengths), \
            "All oligo lengths should be within the specified range"

    def test_handles_short_reference(self):
        reference = "ATGC"
        sequences = ["ATGC", "ATGC"]
        enzyme = Enzyme(fwd_recognition_site="GAATTC", OH_length=4, spacer_length=2)
        min_oligo_len = 2
        max_oligo_len = 25

        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # For short references, only one split is possible
        assert result == [0, len(reference)], "For short references, split indices should be [0, len(reference)]"

    def test_respects_overhangs(self, test_data):
        reference, sequences, enzyme, min_oligo_len, max_oligo_len = test_data

        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # Check that splits align with valid overhang constraints
        valid_overhangs = find_valid_overhangs(reference, range(len(reference)-enzyme.OH_length), enzyme.OH_length)
        assert all(idx in valid_overhangs or idx == 0 or idx == len(reference) for idx in result), \
            "All split indices should align with valid overhangs"


class TestFindValidOverhangs:
    """Integration test suite for the `find_valid_overhangs` function."""

    def test_valid_overhangs_basic(self):
        """Test with a basic valid reference and valid overhang indices."""
        reference = "ATCGGAAACGCGTACGGT"
        OH_compatible_indices = [2, 5, 7]
        OH_len = 3
        expected_result = {
            7: "ACG",  # Reference slice from index 8 to 11
        }

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_all_valid_overhangs(self):
        """Test when no overhangs are valid."""
        reference = "ATCGGATACGCGTACGGT"
        OH_compatible_indices = [1, 4, 7]
        OH_len = 3

        # Assuming none of the overhangs are valid (this will depend on the implementation of `is_valid_overhang`)
        expected_result = {1: 'TCG', 4: 'GAT', 7: 'ACG'}

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_empty_reference(self):
        """Test with an empty reference."""
        reference = ""
        OH_compatible_indices = [0, 1, 2]
        OH_len = 3

        try:
            find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        except ValueError as e:
            assert str(
                e) == "The reference sequence is empty.", f"Expected error message 'The reference sequence is empty.', but got '{str(e)}'"

    def test_no_compatible_indices(self):
        """Test when there are no compatible indices."""
        reference = "ATCGGATACGCGTACGGT"
        OH_compatible_indices = []
        OH_len = 3

        # No overhangs should be found
        expected_result = {}

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_edge_case_overhang(self):
        """Test when the overhang extends beyond the reference sequence."""
        reference = "ATCG"
        OH_compatible_indices = [0, 1, 2]
        OH_len = 5  # Longer than the reference length

        # Expect an IndexError due to out-of-bounds indexing
        with pytest.raises(IndexError, match=r"Index .* with overhang length .* exceeds reference bounds"):
            find_valid_overhangs(reference, OH_compatible_indices, OH_len)

    def test_edge_case_empty_compatible_indices(self):
        """Test when OH_compatible_indices contains empty indices."""
        reference = "ATCGGATACGCGTACGGT"
        OH_compatible_indices = [10, 15]  # Indices that are at the end of the reference
        OH_len = 3

        # Since the indices are near the end, some overhangs might be invalid
        expected_result = {
            10: "CGT",
            15: "GGT"# Valid slice (from 10 to 13)
        }

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"


class TestFindValidCombinations:
    """Integration tests for the `find_valid_combinations` function."""

    def test_basic_case(self):
        """Test a basic case with valid parameters."""
        reference = "AATCGAAAGCCGTGAGTGA"
        OH_compatible_indices = [3, 6, 9, 12]
        min_splits = 1
        max_splits = 4
        true_min_oligo_len = 4
        true_max_oligo_len = 15
        OH_len = 3

        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len)
        # Expect the first valid combination that satisfies all conditions
        expected_result = [9]  # Change to reflect the first valid combination
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_no_valid_combinations(self):
        """Test when no valid combinations exist."""
        reference = "ATGCGTACGTAGCTAGC"
        OH_compatible_indices = [3, 6, 9]
        min_splits = 2
        max_splits = 3
        true_min_oligo_len = 10  # Too large for any valid combination
        true_max_oligo_len = 12
        OH_len = 1

        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len)
        assert result == [], f"Expected an empty list, but got {result}"

    def test_reference_too_short(self):
        """Test when the reference is too short for any splits."""
        reference = "ATGC"
        OH_compatible_indices = [1, 2]
        min_splits = 2
        max_splits = 3
        true_min_oligo_len = 3
        true_max_oligo_len = 4
        OH_len = 1

        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len)
        assert result == [], f"Expected an empty list, but got {result}"

    def test_invalid_indices(self):
        """Test when OH_compatible_indices include invalid positions."""
        reference = "ATGCGTACGTAGCTAGC"
        OH_compatible_indices = [3, 20]  # Index 20 is out of bounds
        min_splits = 1
        max_splits = 2
        true_min_oligo_len = 4
        true_max_oligo_len = 8
        OH_len = 1

        with pytest.raises(ValueError,
                           match="One or more split indices are out of the range of the reference sequence length."):
            find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                    true_max_oligo_len, OH_len)

    def test_empty_reference(self):
        """Test when the reference sequence is empty."""
        reference = ""
        OH_compatible_indices = [0, 1]
        min_splits = 1
        max_splits = 2
        true_min_oligo_len = 4
        true_max_oligo_len = 8
        OH_len = 1

        with pytest.raises(ValueError,
                           match="One or more split indices are out of the range of the reference sequence length."):
            find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                    true_max_oligo_len, OH_len)

    def test_empty_indices(self):
        """Test when OH_compatible_indices is empty."""
        reference = "ATGCGTACGTAGCTAGC"
        OH_compatible_indices = []
        min_splits = 1
        max_splits = 2
        true_min_oligo_len = 4
        true_max_oligo_len = 8
        OH_len = 1

        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                         true_max_oligo_len, OH_len)
        assert result == [], f"Expected an empty list, but got {result}"


class TestCalculateHammingDistances:
    """Test class for calculate_hamming_distances function."""

    def test_no_overlaps(self):
        reference = "ATCCGT"
        OH_len = 3
        split_indices = [0, 3]

        result = calculate_hamming_distances(split_indices, reference, OH_len)
        expected = [3]  # Hamming distances between [0,4], [0,8], [4,8]
        assert result == expected, f"Expected {expected}, but got {result}"

    def test_single_overlap(self):
        """Test when there is a single overlap, the Hamming distance should be computed."""
        reference = "ATCAGT"
        OH_len = 3
        split_indices = [0, 3]

        result = calculate_hamming_distances(split_indices, reference, OH_len)
        expected = [2]  # Hamming distances between [0,3], [3,6]
        assert result == expected, f"Expected {expected}, but got {result}"

    def test_out_of_bounds_index(self):
        """Test when an index is out of bounds, a ValueError should be raised."""
        reference = "ATCGATCGATCG"
        OH_len = 3
        split_indices = [0, 3, 15]  # Index 15 is out of bounds

        with pytest.raises(ValueError, match="One or more indices are out of range of the reference sequence length."):
            calculate_hamming_distances(split_indices, reference, OH_len)

    def test_empty_split_indices(self):
        """Test when there are no split indices, the result should be an empty list."""
        reference = "ATCGATCGATCG"
        OH_len = 2
        split_indices = []

        result = calculate_hamming_distances(split_indices, reference, OH_len)
        expected = []
        assert result == expected, f"Expected {expected}, but got {result}"

    def test_single_split_index(self):
        """Test when there is only a single split index, the result should be an empty list."""
        reference = "ATCGATCGATCG"
        OH_len = 2
        split_indices = [0]

        result = calculate_hamming_distances(split_indices, reference, OH_len)
        expected = []
        assert result == expected, f"Expected {expected}, but got {result}"

    def test_identical_overhangs(self):
        """Test when all overhangs are identical, the Hamming distance should be 0."""
        reference = "ATGATG"
        OH_len = 3
        split_indices = [0, 3]

        result = calculate_hamming_distances(split_indices, reference, OH_len)
        expected = [0]  # Hamming distances between identical overhangs
        assert result == expected, f"Expected {expected}, but got {result}"

    def test_three_different_overhangs(self):
        """Test when the overhangs are different, the Hamming distance should be greater than 0."""
        reference = "ATGGCTATT"
        OH_len = 3
        split_indices = [0, 3, 6]

        result = calculate_hamming_distances(split_indices, reference, OH_len)
        expected = [3, 1, 2]  # Hamming distances between [0,4], [0,8], [4,8]
        assert result == expected, f"Expected {expected}, but got {result}"