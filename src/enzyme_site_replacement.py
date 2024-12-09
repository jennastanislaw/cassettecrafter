"""
Script for identifying and removing enzyme recognition sites from DNA sequences.

This script provides functions to:
1. Remove recognition sites for a specified enzyme from DNA sequences in a DataFrame.
2. Replace individual recognition sites in a DNA sequence using synonymous codons while maintaining the sequence's
functionality.

Modules:
    enzyme_site_replacement_utils: Utility functions for enzyme site matching and codon replacement.

Functions:
    remove_enzyme_sites(enzyme, DNA_df):
        Removes enzyme recognition sites from DNA sequences in a DataFrame.

    replace_enzyme_site(enzyme, DNA):
        Replaces recognition sites for the specified enzyme in a DNA sequence.
"""
import enzyme_site_replacement_utils


def remove_enzyme_sites(enzyme, DNA_df):
    """
    Removes recognition sites for the specified enzyme from DNA sequences in the DataFrame.

    Args:
        enzyme (str): The enzyme whose recognition sites should be removed.
        DNA_df (pd.DataFrame): DataFrame with a 'valid_mixed_bases' column containing DNA sequences.

    Returns:
        pd.DataFrame: Updated DataFrame with a new 'rec_sites_removed' column.
    """
    # Apply the replace_enzyme_site function to each sequence and add results to a new column.
    DNA_df['rec_sites_removed'] = DNA_df['valid_mixed_bases'].apply(lambda dna: replace_enzyme_site(enzyme, dna))

    return DNA_df


def replace_enzyme_site(enzyme, DNA):
    """
    Replaces recognition sites for the specified enzyme in a DNA sequence using synonymous codons.

    Args:
        enzyme (object): Enzyme object containing recognition site details.
        DNA (str): DNA sequence where recognition sites should be replaced.

    Returns:
        str: Modified DNA sequence with recognition sites removed.
    """
    # Convert the DNA sequence to a mutable list for in-place modifications.
    DNA_rec_sites_rem = list(DNA)

    # Find all forward and reverse recognition site matches in the DNA sequence.
    fwd_matches, rev_matches = enzyme_site_replacement_utils.find_matching_sites(enzyme, DNA)

    # Determine codons affected by forward and reverse recognition sites.
    fwd_codons = enzyme_site_replacement_utils.get_affected_codons_by_recognition_sites(
        DNA, fwd_matches, enzyme.fwd_recognition_site)
    rev_codons = enzyme_site_replacement_utils.get_affected_codons_by_recognition_sites(
        DNA, rev_matches, enzyme.rev_recognition_site)

    # Combine all affected codons from forward and reverse matches into a single dictionary.
    all_codons = {**fwd_codons, **rev_codons}

    # Iterate over all recognition sites to replace them with synonymous codons.
    for key, value in all_codons.items():
        # Extract indices of codons associated with the recognition site.
        index_vals = [codon_tuple[1] for codon_tuple in value]
        flattened_index_values = [index for sublist in index_vals for index in sublist]
        max_idx = max(flattened_index_values)
        min_idx = min(flattened_index_values)
        replaced = False  # Track whether a replacement was successfully made.

        # Checks each codon for a synonymous replacement that removes the recognition site
        # As soon as it finds a suitable replacement, exits loop and returns edited sequence
        for codon_tuple in value:
            codon = codon_tuple[0]  # Original codon sequence.
            codon_idx = codon_tuple[1][0]  # Starting index of the codon in the DNA sequence.
            # Generate synonymous codons for the current codon.
            alt_codons = enzyme_site_replacement_utils.generate_synonymous_codons_dna(codon)

            for alt_codon in alt_codons:
                # Replace the current codon with a synonymous codon.
                DNA_rec_sites_rem[codon_idx:codon_idx+3] = alt_codon

                # Check if the recognition site was successfully removed.
                DNA_check = ''.join(DNA_rec_sites_rem)
                fwd_matches, rev_matches = enzyme_site_replacement_utils.find_matching_sites(enzyme, DNA_rec_sites_rem)
                match_idx = fwd_matches + rev_matches
                if key not in match_idx:
                    # Replacement successful; log and exit the loop.
                    print(f"Found replacement for codon {codon} at recognition site index {key}, replaced with {alt_codon}")
                    replaced = True
                    break

            if replaced:  # Stop further replacement attempts for this site.
                break

        # If no suitable replacement was found, raise an error.
        if not replaced:
            raise ValueError(f"No suitable replacement found for recognition site at index {key}.")

    # Join the modified DNA sequence back into a string.
    return ''.join(DNA_rec_sites_rem)

