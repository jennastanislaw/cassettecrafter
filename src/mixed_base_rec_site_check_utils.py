
from itertools import product

from src.enzyme_site_replacement_utils import find_matching_sites

# A dictionary mapping IUPAC codes to possible base sets
IUPAC_to_bases = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],  # Purines
    'Y': ['C', 'T'],  # Pyrimidines
    'S': ['G', 'C'],  # Strong (G or C)
    'W': ['A', 'T'],  # Weak (A or T)
    'K': ['G', 'T'],  # Keto (G or T)
    'M': ['A', 'C'],  # Amino (A or C)
    'B': ['C', 'G', 'T'],  # Not A
    'D': ['A', 'G', 'T'],  # Not C
    'H': ['A', 'C', 'T'],  # Not G
    'V': ['A', 'C', 'G'],  # Not T
    'N': ['A', 'C', 'G', 'T'],  # Any base
}

def expand_dna_sequence(dna_sequence):
    """
    Generate all possible sequences from a DNA sequence containing IUPAC codes.

    Args:
        dna_sequence (str): The DNA sequence containing IUPAC codes.

    Returns:
        list: A list of all possible DNA sequences as strings.
    """
    # Convert each base in the DNA sequence to its possible bases
    base_combinations = []
    for base in dna_sequence:
        if base in IUPAC_to_bases:
            base_combinations.append(IUPAC_to_bases[base])
        else:
            raise ValueError(f"Invalid base '{base}' in DNA sequence.")

    # Generate all possible combinations using itertools.product
    expanded_sequences = [''.join(combination) for combination in product(*base_combinations)]

    return expanded_sequences

def find_non_canonical_bases(sequence):
    """
    Identify the indices of non-canonical bases in a DNA sequence.

    Args:
        sequence (str): A DNA sequence.

    Returns:
        list: A list of tuples where each tuple contains the index and the non-canonical base.
    """
    non_canonical_indices = [
        (idx, base) for idx, base in enumerate(sequence) if base.upper() not in "ACGT"
    ]
    return non_canonical_indices

def check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length):
    """
    Checks if any non-canonical base indices overlap with enzyme recognition sites.

    Args:
        indices (list of tuple): List of (index, base) for non-canonical bases.
        fwd_sites (list of int): Starting indices of forward recognition sites.
        rev_sites (list of int): Starting indices of reverse recognition sites.
        enzyme_oh_length (int): Length of the enzyme's overhang.

    Returns:
        bool: True if any non-canonical base contributes to a recognition site, False otherwise.
    """
    for idx, base in indices:
        for site in fwd_sites + rev_sites:  # Combine forward and reverse sites
            site_range = range(site, site + enzyme_oh_length)

            if idx in site_range:
                print(f"Non-canonical base '{base}' at index {idx} contributes to a recognition site at {site}.")
                return True  # Return True immediately upon finding an overlap

    return False  # Return False if no overlap is found

def check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme):
    """
    Check if any non-canonical base in the expanded sequences overlaps with enzyme recognition sites.
    """
    test = enzyme.fwd_recognition_site

    rec_sites = []
    for seq in expanded_sequences:
        fwd, rev = find_matching_sites(enzyme, seq)
        rec_sites.append(check_non_canonical_in_rec_sites(indices, fwd, rev, len(enzyme.fwd_recognition_site)))
    return rec_sites

def append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows):
    """
    Append sequences with valid recognition sites and non-canonical base info to new rows.
    """
    if any(rec_sites):
        for expanded_seq in expanded_sequences:
            new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': expanded_seq})
    else:
        new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': seq})
