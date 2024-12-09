"""
This script processes DNA sequences within a DataFrame by expanding sequences with IUPAC codes and
validating whether the expanded sequences contain enzyme recognition sites. Valid sequences are appended as new
rows with the original name.

Dependencies:
    pandas (pd): Required for handling DataFrames and managing sequence data.
    mixed_base_rec_site_check_utils (expand_dna_sequence, check_recognition_sites_in_expanded_sequences,
                                      append_valid_sequences, find_IUPAC_codes): Helper functions for
                                      handling DNA sequence expansion and recognition site checks.

Functions:
    degen_codon_checker(df, enzyme):
        Checks and expands DNA sequences with IUPAC codes in the provided DataFrame.
        Appends valid sequences (original and expanded) to a new DataFrame, retaining the original sequence's name.
"""


import pandas as pd
from mixed_base_rec_site_check_utils import (
    expand_dna_sequence,
    check_recognition_sites_in_expanded_sequences,
    append_valid_sequences,
    find_IUPAC_codes
)


def degen_codon_checker(df, enzyme):
    """
    Validate and expand DNA sequences with IUPAC codes in a DataFrame.
    Append expanded sequences as new rows with the original name.

    Parameters:
    df (pd.DataFrame): Input DataFrame with columns 'name' and 'DNA'.
    enzyme (Enzyme): Enzyme object with recognition site properties.

    Returns:
    updated_df (pd.DataFrame): Updated DataFrame with each sequence added as a new row.
    """
    # Create a list to collect new rows
    new_rows = []
    test = enzyme.fwd_recognition_site
    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        seq = row['DNA']
        name = row['name']

        # Check for IUPAC codes
        if any(base not in "ACGT" for base in seq):
            # Find indices of IUPAC codes
            indices = find_IUPAC_codes(seq)

            # Expand sequences based on IUPAC codes
            expanded_sequence_list = expand_dna_sequence(seq)

            # Check recognition sites in expanded sequences
            rec_sites = check_recognition_sites_in_expanded_sequences(indices, expanded_sequence_list, enzyme)

            # Append valid sequences to new rows - if IUPAC code can create a recognition site, the expanded
            # list of sequences with that base are added, otherwise the original sequence containing the IUPAC code
            # is appended
            append_valid_sequences(name, seq, expanded_sequence_list, rec_sites, new_rows)

        else:
            # If no IUPAC codes, add the original sequence as a valid sequence
            new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': seq})

    # Create a DataFrame from the collected new rows
    updated_df = pd.DataFrame(new_rows)

    return updated_df
