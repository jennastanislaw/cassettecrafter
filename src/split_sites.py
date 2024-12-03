from split_sites_utils import find_constant_indices
from split_sites_utils import find_starts_of_consecutive_indices
import pandas as pd
from itertools import combinations

from enzyme_site_replacement_utils import load_enzymes_from_csv
from enzyme_site_replacement_utils import create_enzyme_dict

def generate_cassettes(DNA_df, split_sites):

    for i in range(len(split_sites)-1):

        start = split_sites[i]
        end = split_sites[i + 1]

        DNA_df[f'Cassette {i+1}'] = DNA_df['rec_sites_removed'].apply(lambda x: x[start:end])

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

def find_split_indices(reference, sequences, min_oligo_len, max_oligo_len, enzyme):
    """
    Finds split indices in the reference sequence based on constraints.

    Parameters:
        reference (str): Reference sequence.
        sequences (list): List of aligned sequences.
        min_oligo_len (int): Minimum oligo length.
        max_oligo_len (int): Maximum oligo length.

    Returns:
        pd.DataFrame: DataFrame of split indices, overhangs, and Hamming distances.
    """
    # Find conserved indices
    conserved_indices = find_constant_indices(reference, sequences)

    # Define overhang length (OH_len)
    OH_len = enzyme.OH_length  # Placeholder, replace with value from the sheet if needed

    # Find indices compatible with overhang length
    OH_compatible_indices = find_starts_of_consecutive_indices(conserved_indices, OH_len)

    # Filter indices with valid overhangs based on rules
    valid_overhangs = {
        idx: reference[idx:idx + OH_len]
        for idx in OH_compatible_indices
        if is_valid_overhang(reference[idx:idx + OH_len])
    }

    # Enzyme properties (placeholders; update based on sheet)
    recognition_site_len = len(enzyme.fwd_recognition_site)
    spacer_bps = enzyme.spacer_length

    # Adjust oligo length range
    len_adj = 2 * (recognition_site_len + spacer_bps + OH_len)
    true_min_oligo_len = max(1, min_oligo_len - len_adj)  # Ensure non-negative
    true_max_oligo_len = max_oligo_len - len_adj

    # Determine the range of split site counts
    max_splits = len(reference) // true_min_oligo_len
    min_splits = len(reference) // true_max_oligo_len

    # Find all valid combinations of split indices
    valid_combinations = []
    split = []

    for r in range(min_splits, max_splits):  # Allow splits from 1 to maximum possible
        for split_indices in combinations(OH_compatible_indices, r):
            split_indices = sorted(split_indices)
            oligo_lengths = []

            # Calculate distances from start to each index, and from last index to end
            if split_indices:
                oligo_lengths.append(split_indices[0])  # Distance from start to first index
                oligo_lengths.extend(split_indices[i + 1] - split_indices[i] for i in range(len(split_indices) - 1))
                oligo_lengths.append(len(reference) - split_indices[-1])  # Distance from last index to end

            if all(true_min_oligo_len <= length <= true_max_oligo_len for length in oligo_lengths):
                valid_combinations.append(split_indices)
                valid_combinations.append([0])
                valid_combinations.append([len(reference) - OH_len])

                # Now check the Hamming distances between all pairs of overhangs in this valid combination
                hamming_distances = []

                # Include the first and last overhangs in the split indices
                extended_split_indices = [0] + split_indices + [len(reference) - OH_len]

                for i in range(len(extended_split_indices) - 1):
                    for j in range(i + 1, len(extended_split_indices)):
                        # Get the overhang for the current split index and the next split index
                        overhang1 = reference[extended_split_indices[i]:extended_split_indices[i] + OH_len]
                        overhang2 = reference[extended_split_indices[j]:extended_split_indices[j] + OH_len]

                        # Calculate the Hamming distance between the two overhangs
                        hamming_dist = hamming_distance(overhang1, overhang2)
                        hamming_distances.append(hamming_dist)

                if min(hamming_distances) > 1: #ideally this would be 2 but its actually pretty hard to find overhangs that work with this
                    split = split_indices
                    break

                else:
                    valid_combinations = []

                    # Store the result for each valid combination
                    print(f"Valid Combination: {split_indices}")
                    print(f"Hamming Distances: {hamming_distances}")

            if split:
                break
        if split:
            break
    # No need to return a DataFrame here unless required

    split = [0] + split + [len(reference)]
    return split  # or return a detailed result if necessary

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
