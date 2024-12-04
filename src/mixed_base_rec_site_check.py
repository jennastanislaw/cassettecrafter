from mixed_base_rec_site_check_utils import expand_dna_sequence, check_recognition_sites_in_expanded_sequences, \
    append_valid_sequences
from mixed_base_rec_site_check_utils import find_non_canonical_bases


import pandas as pd

def degen_codon_checker(df, enzyme):
    """
    Validate and expand DNA sequences with non-canonical bases in a DataFrame.
    Append expanded sequences as new rows with the original name.

    Parameters:
    df (pd.DataFrame): Input DataFrame with columns 'name' and 'DNA'.
    enzyme (Enzyme): Enzyme object with recognition site properties.

    Returns:
    pd.DataFrame: Updated DataFrame with each sequence added as a new row.
    """
    # Create a list to collect new rows
    new_rows = []
    test = enzyme.fwd_recognition_site
    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        seq = row['DNA']
        name = row['name']

        # Check for non-standard bases
        if any(base not in "ACGT" for base in seq):
            # Find indices of non-canonical bases
            indices = find_non_canonical_bases(seq)

            # Expand sequences based on non-canonical bases
            expanded_sequence_list = expand_dna_sequence(seq)

            # Check recognition sites in expanded sequences
            rec_sites = check_recognition_sites_in_expanded_sequences(indices, expanded_sequence_list, enzyme)

            # Append valid sequences to new rows
            append_valid_sequences(name, seq, expanded_sequence_list, rec_sites, new_rows)

        else:
            # If no non-standard bases, add the original sequence as a valid sequence
            new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': seq})

    # Create a DataFrame from the collected new rows
    updated_df = pd.DataFrame(new_rows)

    return updated_df
