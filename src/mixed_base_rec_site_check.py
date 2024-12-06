from mixed_base_rec_site_check_utils import expand_dna_sequence
from mixed_base_rec_site_check_utils import find_non_canonical_bases
from mixed_base_rec_site_check_utils import check_non_canonical_in_rec_sites


from enzyme_site_replacement_utils import find_matching_sites
from enzyme_site_replacement_utils import load_enzymes_from_csv
from enzyme_site_replacement_utils import create_enzyme_dict

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

            # Check for recognition sites in expanded sequences
            rec_sites = []
            for seq2 in expanded_sequence_list:
                fwd, rev = find_matching_sites(enzyme, seq2)
                rec_sites.append(check_non_canonical_in_rec_sites(indices, fwd, rev, enzyme.OH_length))

            # If any recognition site issues, append expanded sequences as new rows
            if any(rec_sites):
                for expanded_seq in expanded_sequence_list:
                    new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': expanded_seq})
            else:
                # If no non-standard bases, add the original sequence as a valid sequence
                new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': seq})
        else:
            # If no non-standard bases, add the original sequence as a valid sequence
            new_rows.append({'name': name, 'DNA': seq, 'valid_mixed_bases': seq})

    # Create a DataFrame from the collected new rows
    updated_df = pd.DataFrame(new_rows)

    return updated_df

# Example DataFrame
data = {
    'name': ['seq1', 'seq2', 'seq3'],
    'DNA': ['ATGCN', 'ATGCGTAC', 'ACGTAGAAGRCGCTAG']
}

df = pd.DataFrame(data)

enzyme_class_objs = load_enzymes_from_csv('/Users/siobhan/PycharmProjects/cassettecrafter/src/data/enzyme_sites.csv')

enzyme_dict = create_enzyme_dict(enzyme_class_objs)

enzyme = enzyme_dict.get('BbsI')

print(degen_codon_checker(df, enzyme))
