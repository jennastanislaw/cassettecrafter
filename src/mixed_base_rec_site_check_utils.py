"""
This script processes DNA sequences with IUPAC codes by expanding ambiguous bases into all possible combinations
and checking whether the resulting sequences contain enzyme recognition sites. Sequences with IUPAC codes that do not
create recognition sites are appended to a new list in their original forms, with their associated metadata. If a
sequence has an IUPAC code that creates a recognition site, the list of all possible combinations will be appended to
the new list instead.

Dependencies:
    - itertools.product: Used for generating all combinations of bases for sequences with IUPAC codes.
    - src.enzyme_site_replacement_utils (find_matching_sites): Helper function to find enzyme recognition sites
      within sequences.

Functions:
    expand_dna_sequence(dna_sequence):
        Expands a DNA sequence containing IUPAC codes into all possible sequences.

    find_IUPAC_codes(sequence):
        Identifies indices of IUPAC codes in a given DNA sequence.

    check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length):
        Checks if IUPAC code indices overlap with enzyme recognition sites.

    check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme):
        Verifies if expanded DNA sequences contain enzyme recognition sites influenced by IUPAC codes.

    append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows):
        Appends sequences with valid recognition sites to the results list, retaining metadata.
"""

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
        expanded_sequences (list): A list of all possible DNA sequences as strings.
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

def find_IUPAC_codes(sequence):
    """
    Identify the indices of IUPAC codes in a DNA sequence.

    Args:
        sequence (str): A DNA sequence.

    Returns:
        list: A list of tuples where each tuple contains the index and the IUPAC code.
    """
    IUPAC_code_indices = [
        (idx, base) for idx, base in enumerate(sequence) if base.upper() not in "ACGT"
    ]
    return IUPAC_code_indices

def check_IUPAC_code_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length):
    """
    Checks if any IUPAC code indices overlap with enzyme recognition sites.

    Args:
        indices (list of tuple): List of (index, base) for IUPAC codes.
                                 Each tuple contains the position and the IUPAC code.
        fwd_sites (list of int): Starting indices of forward recognition sites.
                                 These represent positions where the enzyme binds in the forward direction.
        rev_sites (list of int): Starting indices of reverse recognition sites.
                                 These represent positions where the enzyme binds in the reverse direction.
        enzyme_oh_length (int): Length of the enzyme's recognition overhang.

    Returns:
        bool: True if any IUPAC codes contribute to a recognition site, False otherwise.
    """
    # Iterate over each IUPAC code index and base
    for idx, base in indices:
        # Combine forward and reverse recognition site indices into a single list
        for site in fwd_sites + rev_sites:
            # Determine the range of indices covered by the enzyme's recognition site
            site_range = range(site, site + enzyme_oh_length)

            # Check if the IUPAC index falls within the recognition site range
            if idx in site_range:
                print(f"IUPAC code '{base}' at index {idx} contributes to a recognition site at {site}.")
                return True  # Exit immediately upon finding an overlap

    # If no overlaps are found, return False
    return False


def check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme):
    """
    Determine whether IUPAC codes in expanded DNA sequences contribute to enzyme recognition sites.

    Args:
        indices (list of tuple): List of (index, base) for IUPAC codes in the original sequence.
        expanded_sequences (list of str): List of DNA sequences generated by expanding IUPAC codes.
        enzyme (Enzyme): Enzyme object with forward and reverse recognition site properties.

    Returns:
        list of bool: A list indicating whether each expanded sequence creates a recognition site.
    """
    rec_sites = []  # Store recognition site checks for each expanded sequence

    for seq in expanded_sequences:
        # Find forward and reverse recognition site positions in the sequence
        fwd, rev = find_matching_sites(enzyme, seq)

        # Check if IUPAC codes overlap with any recognition site
        rec_sites.append(check_IUPAC_code_in_rec_sites(indices, fwd, rev, len(enzyme.fwd_recognition_site)))

    return rec_sites


def append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows):
    """
    Append sequences to the new rows list based on recognition site validity.

    Args:
        name (str): Name or identifier for the DNA sequence.
        seq (str): Original DNA sequence containing possible IUPAC codes.
        expanded_sequences (list of str): List of DNA sequences generated by expanding IUPAC codes.
        rec_sites (list of bool): List indicating whether each expanded sequence forms a recognition site.
        new_rows (list of dict): List to collect rows for the updated DataFrame.

    Returns:
        None: Modifies the `new_rows` list in place by appending valid sequences.
    """
    # Check if any of the expanded sequences form valid recognition sites
    if any(rec_sites):
        # Append all expanded sequences as valid sequences
        for expanded_seq in expanded_sequences:
            new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': expanded_seq})
    else:
        # If no valid recognition sites, append the original sequence
        new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': seq})

