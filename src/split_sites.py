from src.split_sites_utils import calculate_oligo_lengths, adjust_oligo_lengths, \
    generate_cassettes, is_valid_overhang, hamming_distance, find_constant_indices
from src.split_sites_utils import find_starts_of_consecutive_indices
import pandas as pd
from itertools import combinations

from enzyme_site_replacement_utils import (load_enzymes_from_csv, create_enzyme_dict)


def find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme):
    """
    Finds split indices in the reference sequence based on constraints.

    Parameters:
        reference (str): Reference sequence.
        sequences (list): List of aligned sequences.
        min_oligo_len (int): Minimum oligo length.
        max_oligo_len (int): Maximum oligo length.
        enzyme (Enzyme): Enzyme object with properties like recognition site and overhang length.

    Returns:
        list: List of valid split indices.
    """
    # Step 1: Compute constraints
    conserved_indices = find_constant_indices(reference, sequences)
    OH_len = enzyme.OH_length
    recognition_site_len = len(enzyme.fwd_recognition_site)
    spacer_len = enzyme.spacer_length

    true_min_oligo_len, true_max_oligo_len = adjust_oligo_lengths(
        min_oligo_len, max_oligo_len, recognition_site_len, spacer_len, OH_len
    )

    OH_compatible_indices = find_starts_of_consecutive_indices(conserved_indices, OH_len)

    valid_overhangs = find_valid_overhangs(reference, OH_compatible_indices, OH_len)

    # Step 2: Generate valid split combinations
    max_splits = len(reference) // true_min_oligo_len
    min_splits = len(reference) // true_max_oligo_len

    valid_combination = find_valid_combinations(
        reference, valid_overhangs, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len
    )

    return [0] + valid_combination + [len(reference)]


def find_valid_overhangs(reference, OH_compatible_indices, OH_len):
    """Find valid overhangs based on constraints."""
    if not reference:
        raise ValueError("The reference sequence is empty.")

    reference_len = len(reference)

    for idx in OH_compatible_indices:
        if idx < 0 or idx + OH_len > reference_len or idx > reference_len:
            raise IndexError(
                f"Index {idx} with overhang length {OH_len} exceeds reference bounds (length {reference_len})."
            )

    return {
        idx: reference[idx:idx + OH_len]
        for idx in OH_compatible_indices
        if is_valid_overhang(reference[idx:idx + OH_len])
    }

def find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len):
    """Find the first valid split combination."""
    for r in range(min_splits, max_splits):
        for split_indices in combinations(OH_compatible_indices, r):
            split_indices = sorted(split_indices)
            oligo_lengths = calculate_oligo_lengths(reference, split_indices, OH_len)

            # Check if all oligo lengths are within the valid range
            if all(true_min_oligo_len <= length <= true_max_oligo_len for length in oligo_lengths):
                # Check if Hamming distances for this split are acceptable
                extended_split_indices = [0] + split_indices + [len(reference) - OH_len]
                hamming_distances = calculate_hamming_distances(extended_split_indices, reference, OH_len)

                # If the Hamming distances are valid (i.e., min distance > 1), return the split_indices
                if min(hamming_distances) > 1:
                    return split_indices  # Return the first valid combination

    return []  # Return an empty list if no valid combination is found



def calculate_hamming_distances(split_indices, reference, OH_len):
    """Calculate Hamming distances between overhangs in the split indices."""
    hamming_distances = []

    if any(idx < 0 or idx + OH_len > len(reference) for idx in split_indices):
        raise ValueError("One or more indices are out of range of the reference sequence length.")

    for i in range(len(split_indices) - 1):
        for j in range(i + 1, len(split_indices)):
            overhang1 = reference[split_indices[i]:split_indices[i] + OH_len]
            overhang2 = reference[split_indices[j]:split_indices[j] + OH_len]
            hamming_distances.append(hamming_distance(overhang1, overhang2))
    return hamming_distances


# enzyme_class_objs = load_enzymes_from_csv('/Users/siobhan/PycharmProjects/cassettecrafter/src/data/enzyme_sites.csv')
#
# enzyme_dict = create_enzyme_dict(enzyme_class_objs)
#
# enzyme = enzyme_dict.get('BbsI')
#
# # Example DataFrame
# data = {'rec_sites_removed': ['ATGCGAATTCATCGGGAATTCGCGTGAATTC', 'AAGTCGAATTCGGTTAAGTCGAATTC']}
# DNA_df = pd.DataFrame(data)
#
# # Example split_sites list
# split_sites = [0, 5, 10, 15, 20]
#
# # Generate cassettes
# generate_cassettes(DNA_df, split_sites, enzyme)
#
# print(DNA_df)
#
# # Example usage
# reference = "AGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAA"
# sequences = [
#     "AGAGGTGATGAAGGCAGACAAATCGCTCCAGGGCAAACTGGAAAGACTGCTGATTATAA",
#     "TGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAC",
#     "AGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAA"
# ]



# # Find split indices
# split_sites = find_split_indices(reference, sequences, 10, 40, enzyme)
# print("Split sites with overhangs and Hamming distances:")
# print(split_sites)

# edge cases and other stuff to deal with:
# - what if there are no optimal overhangs
# - test different min and max oligo lengths
# - improve modularity
