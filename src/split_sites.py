from src.split_sites_utils import find_constant_indices
from src.split_sites_utils import find_starts_of_consecutive_indices
import pandas as pd
from itertools import combinations

from enzyme_site_replacement_utils import load_enzymes_from_csv
from enzyme_site_replacement_utils import create_enzyme_dict

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

    valid_combinations = find_valid_combinations(
        reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len
    )

    # Step 3: Check combinations for Hamming distance constraints
    split_indices = select_best_combination(valid_combinations, reference, OH_len)

    return [0] + split_indices + [len(reference)]


def adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len):
    """Adjusts oligo length constraints based on enzyme properties."""
    len_adj = 2 * (recognition_site_len + spacer_bps + OH_len)
    true_min_oligo_len = max(1, min_oligo_len - len_adj)
    true_max_oligo_len = max_oligo_len - len_adj
    return true_min_oligo_len, true_max_oligo_len


def find_valid_overhangs(reference, OH_compatible_indices, OH_len):
    """Find valid overhangs based on constraints."""
    return {
        idx: reference[idx:idx + OH_len]
        for idx in OH_compatible_indices
        if is_valid_overhang(reference[idx:idx + OH_len])
    }


def find_valid_combinations(reference, OH_compatible_indices, min_splits, max_splits, true_min_oligo_len, true_max_oligo_len, OH_len):
    """Find all valid split combinations."""
    valid_combinations = []

    for r in range(min_splits, max_splits):
        for split_indices in combinations(OH_compatible_indices, r):
            split_indices = sorted(split_indices)
            oligo_lengths = calculate_oligo_lengths(reference, split_indices)

            if all(true_min_oligo_len <= length <= true_max_oligo_len for length in oligo_lengths):
                valid_combinations.append(split_indices)

    return valid_combinations


def calculate_oligo_lengths(reference, split_indices):
    """Calculate lengths of oligos based on split indices."""
    oligo_lengths = []
    if split_indices:
        oligo_lengths.append(split_indices[0])  # From start to first index
        oligo_lengths.extend(split_indices[i + 1] - split_indices[i] for i in range(len(split_indices) - 1))
        oligo_lengths.append(len(reference) - split_indices[-1])  # From last index to end
    return oligo_lengths


def select_best_combination(valid_combinations, reference, OH_len):
    """Select the best combination of splits based on Hamming distance constraints."""
    for split_indices in valid_combinations:
        extended_split_indices = [0] + split_indices + [len(reference) - OH_len]
        hamming_distances = calculate_hamming_distances(extended_split_indices, reference, OH_len)

        if min(hamming_distances) > 1:
            return split_indices

    return []  # Return an empty list if no valid combination is found


def calculate_hamming_distances(split_indices, reference, OH_len):
    """Calculate Hamming distances between overhangs in the split indices."""
    hamming_distances = []
    for i in range(len(split_indices) - 1):
        for j in range(i + 1, len(split_indices)):
            overhang1 = reference[split_indices[i]:split_indices[i] + OH_len]
            overhang2 = reference[split_indices[j]:split_indices[j] + OH_len]
            hamming_distances.append(hamming_distance(overhang1, overhang2))
    return hamming_distances


def generate_cassettes(DNA_df, split_sites):

    for i in range(len(split_sites)-1):

        start = split_sites[i]
        end = split_sites[i + 1]

        DNA_df[f'Cassette {i+1}'] = DNA_df['rec_sites_removed'].apply(lambda x: x[start:end])

    return DNA_df

# Example DataFrame
data = {'rec_sites_removed': ['ATGCGAATTCATCGGGAATTCGCGTGAATTC', 'AAGTCGAATTCGGTTAAGTCGAATTC']}
DNA_df = pd.DataFrame(data)

# Example split_sites list
split_sites = [0, 5, 10, 15, 20]

# Generate cassettes
generate_cassettes(DNA_df, split_sites)

print(DNA_df)

def is_valid_overhang(overhang):
    """Checks if an overhang meets the defined rules."""
    # Rule 2: Avoid palindromes
    if overhang == overhang[::-1]:
        return False

    # Rule 3: No overhangs with the same three nucleotides in a row
    if any(overhang[i] == overhang[i+1] == overhang[i+2] for i in range(len(overhang) - 2)):
        return False

    # Rule 5: Avoid overhangs with extreme GC content (0% or 100%)
    gc_content = (overhang.count('G') + overhang.count('C')) / len(overhang)
    if gc_content == 0 or gc_content == 1:
        return False

    return True

def hamming_distance(seq1, seq2):
    """Calculates the Hamming distance between two sequences."""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

# Example usage
reference = "AGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAA"
sequences = [
    "AGAGGTGATGAAGGCAGACAAATCGCTCCAGGGCAAACTGGAAAGACTGCTGATTATAA",
    "TGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAC",
    "AGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAA"
]


enzyme_class_objs = load_enzymes_from_csv('/Users/siobhan/PycharmProjects/cassettecrafter/src/data/enzyme_sites.csv')

enzyme_dict = create_enzyme_dict(enzyme_class_objs)

enzyme = enzyme_dict.get('BbsI')
# Find split indices
split_sites = find_split_indices(reference, sequences, 10, 40, enzyme)
print("Split sites with overhangs and Hamming distances:")
print(split_sites)

# edge cases and other stuff to deal with:
# - what if there are no optimal overhangs
# - test different min and max oligo lengths
# - improve modularity
