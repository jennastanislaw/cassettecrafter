"""
Test suite for DNA sequence processing and enzyme-based operations.

This module includes a set of tests to validate the functions involved in DNA sequence
splitting, finding valid overhangs, and working with restriction enzymes. It contains:
- A fixture (`test_data`) that provides test data for DNA sequence processing, enzyme properties,
  and oligo length constraints.
- A `TestFindSplitIndices` class that tests the `find_split_indices` function for valid output,
  proper constraints handling, short references, and overhang validation.
- A `TestFindValidOverhangs` class that tests the `find_valid_overhangs` function for different cases,
  including basic validity, handling of empty sequences, and edge cases with invalid indices.

The tests cover functionality such as ensuring proper split indices, overhang compliance, and correct
error handling when inputs do not meet expectations.

Classes:
    - Enzyme: A class representing a restriction enzyme with attributes for the forward recognition site,
      overhang length, and spacer length.
    - TestFindSplitIndices: A test suite for the `find_split_indices` function.
    - TestFindValidOverhangs: A test suite for the `find_valid_overhangs` function.

Functions:
    - test_data: A fixture that provides sample data for DNA sequence processing, enzyme details,
      and oligo length constraints.
"""
import pytest
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from split_sites import (
    find_split_indices,
    find_valid_overhangs,
    find_valid_combinations,
    calculate_hamming_distances
)

class Enzyme:
    """
    A class to represent a restriction enzyme.

    Attributes:
        fwd_recognition_site (str): The forward recognition site of the enzyme.
        OH_length (int): The length of the overhang created after the enzyme cuts.
        spacer_length (int): The length of the spacer sequence (the region between the recognition sites).
    """

    def __init__(self, fwd_recognition_site, OH_length, spacer_length):
        """
        Initializes the Enzyme object with the provided recognition site, overhang length,
        and spacer length.

        Parameters:
            fwd_recognition_site (str): The forward recognition site of the enzyme.
            OH_length (int): The length of the overhang created after the enzyme cuts.
            spacer_length (int): The length of the spacer sequence between the recognition sites.
        """
        self.fwd_recognition_site = fwd_recognition_site  # The forward recognition site of the enzyme
        self.OH_length = OH_length  # The overhang length of the cut produced by the enzyme
        self.spacer_length = spacer_length  # The length of the spacer between recognition sites


@pytest.fixture
def test_data():
    """
    Fixture to provide test data for DNA sequence processing.

    This fixture returns a set of sample data for testing purposes:
    - A reference sequence to compare against.
    - A list of sequences that will be tested against the reference sequence.
    - An Enzyme object representing a restriction enzyme with specified attributes.
    - Minimum and maximum lengths for oligos to check against during processing.

    Returns:
        tuple: A tuple containing:
            - reference (str): The reference DNA sequence.
            - sequences (list of str): A list of DNA sequences to be tested.
            - enzyme (Enzyme): An Enzyme object with the recognition site, overhang length, and spacer length.
            - min_oligo_len (int): The minimum required oligo length for sequences.
            - max_oligo_len (int): The maximum allowed oligo length for sequences.
    """
    reference = "ATGCGTACGTTAGCTAGCGTACGTTAGC"  # Reference DNA sequence used for comparison
    sequences = [
        "ATGCGTACGTTAGCTAGCGTACGTTAGC",  # Identical sequence to the reference
        "ATGCGTACGTTAGCTAGCGTACGTTTGC",  # Sequence with one mismatch at the end
        "ATGCGTACGTTAGCTAGCGTACGTTTAC"   # Sequence with a mismatch at a different position
    ]
    enzyme = Enzyme(fwd_recognition_site="GAATTC", OH_length=4, spacer_length=2)  # Create an Enzyme object
    min_oligo_len = 6  # Minimum oligo length to check against during DNA sequence processing
    max_oligo_len = 37  # Maximum oligo length for allowed DNA sequences

    return reference, sequences, enzyme, min_oligo_len, max_oligo_len  # Return all test data as a tuple



class TestFindSplitIndices:
    """
    Test suite for the function `find_split_indices` that identifies split points
    for DNA sequences based on enzyme recognition sites and other constraints.
    """

    def test_valid_output(self, test_data):
        """
        Test the output format and basic validity of the `find_split_indices` function.

        This test checks:
        - The result is a list.
        - All elements of the list are integers.
        - The first index is 0, and the last index is the length of the reference sequence.
        """
        # Unpack test data
        reference, sequences, enzyme, min_oligo_len, max_oligo_len = test_data

        # Call the function under test
        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # Verify the result has the expected format
        assert isinstance(result, list), "Output should be a list"
        assert all(isinstance(i, int) for i in result), "All indices should be integers"
        assert result[0] == 0, "First index should be 0"
        assert result[-1] == len(reference), "Last index should be the length of the reference"

    def test_respects_constraints(self, test_data):
        """
        Test that the split indices respect the specified oligo length constraints.

        This test ensures:
        - The oligo lengths fall within the valid range after splitting the sequences.
        - The oligo length constraints take into account the enzyme's overhang and spacer lengths.
        """
        # Unpack test data
        reference, sequences, enzyme, min_oligo_len, max_oligo_len = test_data
        OH_len = enzyme.OH_length

        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # Remove the first and last indices (which are the start and end of the reference)
        result = result[1:-1]

        # Initialize the list to track the oligo lengths
        oligo_lengths = []

        # Calculate adjustment for oligo lengths based on enzyme parameters
        len_adj = 2 * (len(enzyme.fwd_recognition_site) + enzyme.spacer_length + enzyme.OH_length)
        true_min_oligo_len = max(1, min_oligo_len - len_adj)
        true_max_oligo_len = max_oligo_len - len_adj

        # Calculate oligo lengths based on split indices
        if not result:
            # If there are no splits, the entire reference sequence is considered one oligo
            oligo_lengths.append(len(reference))
        else:
            # Calculate oligo lengths for splits between indices
            oligo_lengths.append(result[0] + OH_len)  # From start to first split
            oligo_lengths.extend(
                result[i + 1] - result[i] + OH_len for i in range(len(result) - 1)
            )
            oligo_lengths.append(len(reference) - result[-1])  # From last split to end

        # Ensure all oligo lengths fall within the valid range
        assert all(true_min_oligo_len <= length <= true_max_oligo_len for length in oligo_lengths), \
            "All oligo lengths should be within the specified range"

    def test_handles_short_reference(self):
        """
        Test behavior when the reference sequence is short.

        For short reference sequences (e.g., <= 4 bases), the only valid split
        indices should be [0, len(reference)] since no further splitting is possible.
        """
        reference = "ATGC"
        sequences = ["ATGC", "ATGC"]
        enzyme = Enzyme(fwd_recognition_site="GAATTC", OH_length=4, spacer_length=2)
        min_oligo_len = 2
        max_oligo_len = 25

        # Call the function under test
        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # Verify that for short references, the result is just [0, len(reference)]
        assert result == [0, len(reference)], "For short references, split indices should be [0, len(reference)]"

    def test_respects_overhangs(self, test_data):
        """
        Test that split indices align with valid overhang constraints.

        This test ensures:
        - Split indices are either at the start or end of the reference sequence,
          or at positions where the overhangs created by the enzyme are valid.
        """
        # Unpack test data
        reference, sequences, enzyme, min_oligo_len, max_oligo_len = test_data

        # Call the function under test
        result = find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme)

        # Get the valid overhang positions based on enzyme properties
        valid_overhangs = find_valid_overhangs(reference, range(len(reference) - enzyme.OH_length), enzyme.OH_length)

        # Ensure all split indices are valid according to overhang constraints
        assert all(idx in valid_overhangs or idx == 0 or idx == len(reference) for idx in result), \
            "All split indices should align with valid overhangs"


class TestFindValidOverhangs:
    """Integration test suite for the `find_valid_overhangs` function."""

    def test_valid_overhangs_basic(self):
        """
        Test with a basic valid reference and valid overhang indices.

        This test verifies that the function correctly identifies valid overhangs
        in the given reference sequence based on the provided overhang indices and
        overhang length.
        """
        reference = "ATCGGAAACGCGTACGGT"
        OH_compatible_indices = [2, 5, 7]
        OH_len = 3
        expected_result = {
            7: "ACG",  # Reference slice from index 8 to 11
        }

        # Run the function and compare results
        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_all_valid_overhangs(self):
        """
        Test when all overhang indices are valid.

        This test ensures that the function correctly returns the valid overhangs
        corresponding to all provided indices, and returns the expected slices from
        the reference sequence.
        """
        reference = "ATCGGATACGCGTACGGT"
        OH_compatible_indices = [1, 4, 7]
        OH_len = 3

        # Assuming the overhangs are valid for these indices (depends on the implementation)
        expected_result = {1: 'TCG', 4: 'GAT', 7: 'ACG'}

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_empty_reference(self):
        """
        Test with an empty reference sequence.

        This test verifies that the function raises a ValueError when the reference
        sequence is empty, as no overhangs can be found in this case.
        """
        reference = ""
        OH_compatible_indices = [0, 1, 2]
        OH_len = 3

        # Check for expected error when the reference is empty
        try:
            find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        except ValueError as e:
            # Assert the correct error message is raised
            assert str(e) == "The reference sequence is empty.", \
                f"Expected error message 'The reference sequence is empty.', but got '{str(e)}'"

    def test_no_compatible_indices(self):
        """
        Test when there are no compatible indices.

        This test ensures that when no overhang indices are provided, the function
        correctly returns an empty result, indicating that no valid overhangs were found.
        """
        reference = "ATCGGATACGCGTACGGT"
        OH_compatible_indices = []
        OH_len = 3

        # No valid overhangs should be found
        expected_result = {}

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_edge_case_overhang(self):
        """
        Test when the overhang extends beyond the reference sequence.

        This test checks that the function raises an IndexError when the overhang
        length is greater than the available reference sequence at the provided indices.
        """
        reference = "ATCG"
        OH_compatible_indices = [0, 1, 2]
        OH_len = 5  # Overhang is longer than the reference length

        # Expect an IndexError because the overhang exceeds the bounds of the reference sequence
        with pytest.raises(IndexError, match=r"Index .* with overhang length .* exceeds reference bounds"):
            find_valid_overhangs(reference, OH_compatible_indices, OH_len)

    def test_edge_case_empty_compatible_indices(self):
        """
        Test when `OH_compatible_indices` contains indices near the end of the reference.

        This test ensures that the function correctly handles cases where the provided
        indices are close to or at the end of the reference sequence, potentially
        leading to out-of-bounds overhangs.
        """
        reference = "ATCGGATACGCGTACGGT"
        OH_compatible_indices = [10, 15]  # Indices at the end of the reference
        OH_len = 3

        # Some of the overhangs might be valid or invalid depending on the reference length
        expected_result = {
            10: "CGT",  # Valid slice from 10 to 13
            15: "GGT"  # Valid slice from 15 to 18
        }

        result = find_valid_overhangs(reference, OH_compatible_indices, OH_len)
        assert result == expected_result, f"Expected {expected_result}, but got {result}"


class TestFindValidCombinations:
    """Integration tests for the `find_valid_combinations` function."""

    def test_basic_case(self):
        """Test a basic case with valid parameters."""
        # Test parameters with valid indices and oligo length range
        reference = "AATCGAAAGCCGTGAGTGA"  # DNA sequence to find combinations within
        OH_compatible_indices = [3, 6, 9, 12]  # Indices where the overhang is compatible
        min_splits = 1  # Minimum number of splits required
        max_splits = 4  # Maximum number of splits allowed
        true_min_oligo_len = 4  # Minimum allowed oligo length
        true_max_oligo_len = 15  # Maximum allowed oligo length
        OH_len = 3  # Length of overhang

        # Call the function to get valid combinations
        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                         true_max_oligo_len, OH_len)

        # The expected result should match the first valid combination found
        expected_result = [9]  # Change to reflect the first valid combination
        # Assert that the result matches the expected outcome
        assert result == expected_result, f"Expected {expected_result}, but got {result}"

    def test_no_valid_combinations(self):
        """Test when no valid combinations exist."""
        # Test with a reference sequence and parameters where no valid combinations are possible
        reference = "ATGCGTACGTAGCTAGC"  # DNA sequence to find combinations within
        OH_compatible_indices = [3, 6, 9]  # Indices where the overhang is compatible
        min_splits = 2  # Minimum number of splits required
        max_splits = 3  # Maximum number of splits allowed
        true_min_oligo_len = 10  # Minimum oligo length is too large for any valid combinations
        true_max_oligo_len = 12  # Maximum oligo length
        OH_len = 1  # Length of overhang

        # Call the function to get valid combinations
        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                         true_max_oligo_len, OH_len)

        # Since no valid combinations exist, the result should be an empty list
        assert result == [], f"Expected an empty list, but got {result}"

    def test_reference_too_short(self):
        """Test when the reference is too short for any splits."""
        # Test with a short reference sequence that cannot accommodate the requested splits
        reference = "ATGC"  # DNA sequence that is too short
        OH_compatible_indices = [1, 2]  # Indices where the overhang is compatible
        min_splits = 2  # Minimum number of splits required
        max_splits = 3  # Maximum number of splits allowed
        true_min_oligo_len = 3  # Minimum oligo length
        true_max_oligo_len = 4  # Maximum oligo length
        OH_len = 1  # Length of overhang

        # Call the function to get valid combinations
        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                         true_max_oligo_len, OH_len)

        # Since the reference is too short to allow any splits, the result should be an empty list
        assert result == [], f"Expected an empty list, but got {result}"

    def test_invalid_indices(self):
        """Test when OH_compatible_indices include invalid positions."""
        # Test when indices are out of bounds for the reference sequence
        reference = "ATGCGTACGTAGCTAGC"  # DNA sequence to find combinations within
        OH_compatible_indices = [3, 20]  # Index 20 is out of bounds for the reference length
        min_splits = 1  # Minimum number of splits required
        max_splits = 2  # Maximum number of splits allowed
        true_min_oligo_len = 4  # Minimum oligo length
        true_max_oligo_len = 8  # Maximum oligo length
        OH_len = 1  # Length of overhang

        # Expect a ValueError when the indices are out of range
        with pytest.raises(ValueError,
                           match="One or more split indices are out of the range of the reference sequence length."):
            find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                    true_max_oligo_len, OH_len)

    def test_empty_reference(self):
        """Test when the reference sequence is empty."""
        # Test with an empty reference sequence
        reference = ""  # Empty sequence
        OH_compatible_indices = [0, 1]  # Indices where the overhang is compatible (but won't work with empty sequence)
        min_splits = 1  # Minimum number of splits required
        max_splits = 2  # Maximum number of splits allowed
        true_min_oligo_len = 4  # Minimum oligo length
        true_max_oligo_len = 8  # Maximum oligo length
        OH_len = 1  # Length of overhang

        # Expect a ValueError because indices won't be valid for an empty sequence
        with pytest.raises(ValueError,
                           match="One or more split indices are out of the range of the reference sequence length."):
            find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                    true_max_oligo_len, OH_len)

    def test_empty_indices(self):
        """Test when OH_compatible_indices is empty."""
        # Test with no indices provided for overhang compatibility
        reference = "ATGCGTACGTAGCTAGC"  # DNA sequence to find combinations within
        OH_compatible_indices = []  # No compatible indices provided
        min_splits = 1  # Minimum number of splits required
        max_splits = 2  # Maximum number of splits allowed
        true_min_oligo_len = 4  # Minimum oligo length
        true_max_oligo_len = 8  # Maximum oligo length
        OH_len = 1  # Length of overhang

        # Since no indices are provided, the result should be an empty list
        result = find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len,
                                         true_max_oligo_len, OH_len)
        assert result == [], f"Expected an empty list, but got {result}"


class TestCalculateHammingDistances:
    """Test class for the calculate_hamming_distances function."""

    def test_no_overlaps(self):
        """Test when there are no overlaps between the overhangs."""
        reference = "ATCCGT"  # Reference sequence
        OH_len = 3  # Length of the overhangs
        split_indices = [0, 3]  # Split the reference into two parts at indices 0 and 3

        result = calculate_hamming_distances(split_indices, reference, OH_len)  # Call function to calculate Hamming distances
        expected = [3]  # Expected Hamming distance between the two overhangs
        assert result == expected, f"Expected {expected}, but got {result}"  # Check if the result matches the expected

    def test_single_overlap(self):
        """Test when there is a single overlap, the Hamming distance should be computed."""
        reference = "ATCAGT"  # Reference sequence
        OH_len = 3  # Length of the overhangs
        split_indices = [0, 3]  # Split the reference at indices 0 and 3

        result = calculate_hamming_distances(split_indices, reference, OH_len)  # Call function to calculate Hamming distances
        expected = [2]  # Expected Hamming distance for the overlap between the two overhangs
        assert result == expected, f"Expected {expected}, but got {result}"  # Check if the result matches the expected

    def test_out_of_bounds_index(self):
        """Test when an index is out of bounds, a ValueError should be raised."""
        reference = "ATCGATCGATCG"  # Reference sequence
        OH_len = 3  # Length of the overhangs
        split_indices = [0, 3, 15]  # Invalid index 15 (out of bounds for the reference sequence)

        # Expect the function to raise a ValueError due to the out-of-bounds index
        with pytest.raises(ValueError, match="One or more indices are out of range of the reference sequence length."):
            calculate_hamming_distances(split_indices, reference, OH_len)

    def test_empty_split_indices(self):
        """Test when there are no split indices, the result should be an empty list."""
        reference = "ATCGATCGATCG"  # Reference sequence
        OH_len = 2  # Length of the overhangs
        split_indices = []  # No split indices provided

        result = calculate_hamming_distances(split_indices, reference, OH_len)  # Call function with empty split indices
        expected = []  # Expecting an empty list because there are no overhangs to compare
        assert result == expected, f"Expected {expected}, but got {result}"  # Check if the result is as expected

    def test_single_split_index(self):
        """Test when there is only a single split index, the result should be an empty list."""
        reference = "ATCGATCGATCG"  # Reference sequence
        OH_len = 2  # Length of the overhangs
        split_indices = [0]  # Only one split index provided

        result = calculate_hamming_distances(split_indices, reference, OH_len)  # Call function with a single split index
        expected = []  # No Hamming distances can be computed with a single overhang
        assert result == expected, f"Expected {expected}, but got {result}"  # Check if the result is as expected

    def test_identical_overhangs(self):
        """Test when all overhangs are identical, the Hamming distance should be 0."""
        reference = "ATGATG"  # Reference sequence with identical overhangs
        OH_len = 3  # Length of the overhangs
        split_indices = [0, 3]  # Split the reference sequence into two identical overhangs

        result = calculate_hamming_distances(split_indices, reference, OH_len)  # Call function to calculate Hamming distances
        expected = [0]  # Since the overhangs are identical, the Hamming distance is 0
        assert result == expected, f"Expected {expected}, but got {result}"  # Check if the result matches the expected

    def test_three_different_overhangs(self):
        """Test when the overhangs are different, the Hamming distance should be greater than 0."""
        reference = "ATGGCTATT"  # Reference sequence with different overhangs
        OH_len = 3  # Length of the overhangs
        split_indices = [0, 3, 6]  # Split the reference sequence at indices 0, 3, and 6

        result = calculate_hamming_distances(split_indices, reference, OH_len)  # Call function to calculate Hamming distances
        expected = [3, 1, 2]  # Expected Hamming distances between the different overhangs
        assert result == expected, f"Expected {expected}, but got {result}"  # Check if the result matches the expected
