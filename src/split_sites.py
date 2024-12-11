"""
This script contains functions for finding valid split indices in a reference DNA sequence
based on enzyme recognition site, oligo length, and overhang sequence constraints. It involves calculating
overhang sequences, checking their validity, and determining appropriate split combinations
while considering Hamming distances between the overhangs.

Functions:
- find_split_indices: Identifies valid split indices for the reference sequence, based on
  enzyme properties, oligo length constraints, and recognition sites.
- find_valid_overhangs: Returns valid overhang sequences based on given indices and
  sequence constraints.
- find_valid_combinations: Finds the first valid combination of split indices that meets
  oligo length and Hamming distance requirements.
- calculate_hamming_distances: Calculates the Hamming distances between all pairwise
  overhangs generated from split indices.

Dependencies:
- `src.split_sites_utils`: Functions for oligo length adjustments, Hamming distance
  calculations, and determining conserved indices.
- `itertools`: Provides combinations function for generating split index combinations.
- `math`: Used for mathematical operations like ceiling of the split index range.

"""

# from src.split_sites_utils import (calculate_oligo_lengths, adjust_oligo_lengths, is_valid_overhang, hamming_distance,
#                                    find_constant_indices)
# from src.split_sites_utils import find_starts_of_consecutive_indices
from utils.split_sites_utils import (calculate_oligo_lengths, adjust_oligo_lengths, is_valid_overhang, hamming_distance, find_starts_of_consecutive_indices, find_constant_indices)
#from split_sites_utils import find_starts_of_consecutive_indices
from itertools import combinations
import math


def find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme):
    """
    Finds valid split indices in the reference sequence based on enzyme recognition and oligo length constraints.

    Args:
        reference (str): The reference DNA sequence.
        sequences (list): List of aligned DNA sequences.
        min_oligo_len (int): Minimum allowed oligo length.
        max_oligo_len (int): Maximum allowed oligo length.
        enzyme (Enzyme): Enzyme object containing recognition site, overhang, and spacer length.

    Returns:
        list: List of indices where the reference sequence can be split.
    """
    # Step 1: Compute constraints
    conserved_indices = find_constant_indices(reference, sequences)  # Find conserved positions across sequences
    OH_len = enzyme.OH_length  # Overhang length
    recognition_site_len = len(enzyme.fwd_recognition_site)  # Recognition site length
    spacer_len = enzyme.spacer_length  # Spacer length

    # Adjust oligo lengths based on enzyme properties
    true_min_oligo_len, true_max_oligo_len = adjust_oligo_lengths(
        min_oligo_len, max_oligo_len, recognition_site_len, spacer_len, OH_len
    )

    # Step 2: Identify overhang-eligible indices
    OH_compatible_indices = find_starts_of_consecutive_indices(conserved_indices, OH_len)
    valid_overhangs = find_valid_overhangs(reference, OH_compatible_indices, OH_len)

    # Step 3: Generate valid split combinations
    max_splits = len(reference) // true_min_oligo_len
    min_splits = math.ceil(len(reference) / true_max_oligo_len)

    # Find valid combinations of splits within the calculated range
    valid_combination = find_valid_combinations(
        reference, valid_overhangs, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len
    )

    # Return valid split indices including the start (0) and end (len(reference))
    return [0] + valid_combination + [len(reference)]



def find_valid_overhangs(reference, OH_compatible_indices, OH_len):
    """
    Generates a dictionary of indices where the overhang that would be generated meets the sequence pattern requirements
    for a golden gate assembly overhang.

    Args:
        reference (str): The reference DNA sequence.
        OH_compatible_indices (list): List of indices where overhangs are possible.
        OH_len (int): The length of the overhang to extract.

    Returns:
        dict: A dictionary with valid overhangs, where keys are indices and values are overhang sequences.

    Raises:
        ValueError: If the reference sequence is empty.
        IndexError: If any overhang index exceeds the bounds of the reference sequence.
    """
    if not reference:
        raise ValueError("The reference sequence is empty.")

    reference_len = len(reference)

    # Check if the indices are within the reference bounds and the overhangs are valid
    for idx in OH_compatible_indices:
        if idx < 0 or idx + OH_len > reference_len or idx > reference_len:
            raise IndexError(
                f"Index {idx} with overhang length {OH_len} exceeds reference bounds (length {reference_len})."
            )

    # Return a dictionary of valid overhangs based on the defined criteria
    return {
        idx: reference[idx:idx + OH_len]
        for idx in OH_compatible_indices
        if is_valid_overhang(reference[idx:idx + OH_len])
    }


def find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len):
    """
    Find the first valid split combination of the reference sequence.

    Args:
        reference (str): The reference DNA sequence.
        OH_compatible_indices (list): List of indices where overhangs are compatible.
        min_splits (int): Minimum number of splits allowed.
        max_splits (int): Maximum number of splits allowed.
        true_min_oligo_len (int): Minimum valid oligo length.
        true_max_oligo_len (int): Maximum valid oligo length.
        OH_len (int): Length of the overhang.

    Returns:
        list: A list of split indices that represent a valid combination.
              If no valid combination is found, an empty list is returned.
    """
    for r in range(min_splits, max_splits):
        # Generate all possible combinations of split indices for the current number of splits
        for split_indices in combinations(OH_compatible_indices, r):
            split_indices = sorted(split_indices)

            # Calculate the lengths of the oligos generated by this split
            oligo_lengths = calculate_oligo_lengths(reference, split_indices, OH_len)

            # Check if all oligo lengths are within the valid range
            if all(true_min_oligo_len <= length <= true_max_oligo_len for length in oligo_lengths):
                # Create extended split indices and calculate Hamming distances
                extended_split_indices = [0] + split_indices + [len(reference) - OH_len]
                hamming_distances = calculate_hamming_distances(extended_split_indices, reference, OH_len)

                # If the Hamming distances are valid (min distance > 1), return this split
                if min(hamming_distances) > 1:
                    return split_indices  # Return the first valid combination

    return []  # Return an empty list if no valid combination is found


def calculate_hamming_distances(split_indices, reference, OH_len):
    """
    Calculate Hamming distances between overhangs for all pairwise splits.

    Args:
        split_indices (list): List of indices where the sequence is split.
        reference (str): The reference DNA sequence.
        OH_len (int): Length of the overhang.

    Returns:
        list: A list of Hamming distances between each pair of overhangs from the split indices.

    Raises:
        ValueError: If any split index or overhang is out of range of the reference sequence.
    """
    hamming_distances = []

    # Check if any split index is out of range of the reference sequence
    if any(idx < 0 or idx + OH_len > len(reference) for idx in split_indices):
        raise ValueError("One or more indices are out of range of the reference sequence length.")

    # Calculate Hamming distances between all pairwise overhangs
    for i in range(len(split_indices) - 1):
        for j in range(i + 1, len(split_indices)):
            overhang1 = reference[split_indices[i]:split_indices[i] + OH_len]
            overhang2 = reference[split_indices[j]:split_indices[j] + OH_len]
            hamming_distances.append(hamming_distance(overhang1, overhang2))

    return hamming_distances

