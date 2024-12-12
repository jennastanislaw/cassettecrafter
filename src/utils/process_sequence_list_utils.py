"""
This script provides functions for manipulating and modifying DNA sequences using enzyme recognition sites.
It includes functionality for generating random DNA sequences, adding enzyme recognition sites with spacers,
ensuring a minimum length for DNA sequences, and checking for the presence of enzyme recognition sites.

Dependencies:
    - Biopython (Bio.Seq): Required for handling and manipulating DNA sequences.
    - NumPy: Used for mathematical operations like calculating sequence lengths.
    - enzyme_site_replacement_utils (custom module): Contains utility functions for enzyme site matching, enzyme
    loading, and creating enzyme dictionaries.

Functions:
    generate_random_dna(n):
        Generates a random DNA sequence of the specified length.

    check_random_oligo_for_sites(DNA, enzyme):
        Checks for the presence of enzyme recognition sites in a given DNA sequence.

    add_end_restriction_sites(df, column, enzyme):
        Adds enzyme recognition sites (with spacers if necessary) to DNA sequences in a DataFrame column.

    generate_spacers_rest_sites(dna, enzyme):
        Adds enzyme recognition sites and spacers to a DNA sequence.

    generate_spacer(dna, rec_site, enzyme, direction="forward"):
        Generates a random spacer and adds it to the DNA sequence either before or after the sequence
        based on the specified direction.

    generate_siteless_sequence(n, enzyme):
        Generates a random DNA sequence that does not contain any enzyme recognition sites.

    ensure_minimum_length(df, col, min_oligo_size, enzyme):
        Ensures that each DNA sequence in a DataFrame column meets the minimum length requirement by
        adding random siteless sequences at both ends.
"""

from Bio.Seq import Seq
import random
import numpy as np
from utils.enzyme_site_replacement_utils import find_matching_sites, load_enzymes_from_csv, create_enzyme_dict


def generate_random_dna(n):
    """
    Generate a random DNA sequence of specified length.

    Args:
        n (int): Length of the DNA sequence to generate. Must be a positive integer.

    Returns:
        Seq: A Biopython Seq object representing the random DNA sequence.

    Raises:
        ValueError: If the specified length `n` is negative.
    """
    if n == 0:
        return Seq("")  # Return an empty DNA sequence for length 0

    if n < 0:
        raise ValueError("Length of DNA sequence must be a positive integer.")

    # Define the DNA nucleotide bases
    nucleotides = ['A', 'T', 'C', 'G']

    # Generate a random sequence of length `n` using the nucleotide bases
    random_sequence = ''.join(random.choices(nucleotides, k=n))

    # Convert the string sequence to a Biopython Seq object
    dna_sequence = Seq(random_sequence)

    return dna_sequence


def check_random_oligo_for_sites(DNA, enzyme):
    """
    Check for the presence of enzyme recognition sites in a given DNA sequence.

    Args:
        DNA (str): DNA sequence to analyze.
        enzyme (Enzyme): Enzyme object containing recognition site properties.

    Returns:
        int: Total number of forward and reverse recognition site matches found in the DNA sequence.
    """
    # Find forward and reverse recognition site matches
    fwd_matches, rev_matches = find_matching_sites(enzyme, DNA)

    # Return the total number of matches
    return len(rev_matches + fwd_matches)


def add_end_restriction_sites(df, column, enzyme):
    """
    Add enzyme recognition sites to DNA sequences in a DataFrame column.

    Prepends the forward recognition site and appends the reverse recognition site
    (with spacers if necessary) to each DNA sequence.

    Args:
        df (pd.DataFrame): DataFrame containing DNA sequences.
        column (str): Name of the column containing DNA sequences to modify.
        enzyme (Enzyme): Enzyme object.

    Returns:
        pd.DataFrame: Updated DataFrame with modified DNA sequences in the specified column.
    """

    # Apply the function to add spacers and recognition sites to each sequence
    df[column] = df[column].apply(lambda dna: generate_spacers_rest_sites(dna, enzyme))

    return df



def generate_spacers_rest_sites(dna, enzyme):
    """
    This function adds a forward enzyme recognition site at the beginning of the DNA sequence,
    followed by a spacer, then the original DNA sequence, followed by a reverse recognition site
    at the end with a spacer.

    Args:
        dna (str): The original DNA sequence.
        enzyme (Enzyme): Enzyme object containing the recognition sites and spacer length.

    Returns:
        str: The modified DNA sequence with recognition sites and spacers added.
    """
    # Extract enzyme recognition sites
    forward_recognition_site = enzyme.fwd_recognition_site
    reverse_recognition_site = enzyme.rev_recognition_site

    # Generate spacers for forward and reverse recognition sites
    forward_rest_spacer = generate_spacer(dna, forward_recognition_site, enzyme, 'forward')
    reverse_rest_spacer = generate_spacer(dna, reverse_recognition_site, enzyme, 'reverse')

    # Construct the new DNA sequence with recognition sites and spacers
    dna_new = forward_recognition_site + forward_rest_spacer + dna + reverse_rest_spacer + reverse_recognition_site

    return str(dna_new)



def generate_spacer(dna, rec_site, enzyme, direction="forward"):
    """
    This function generates a random spacer of a specified length, adds it either before or after
    the DNA sequence (based on the direction), and checks if the resulting sequence meets the
    enzyme recognition site requirements. If the sequence is valid, the spacer is returned.

    Args:
        dna (str): The DNA sequence to which the spacer will be added.
        rec_site (str): The enzyme recognition site (either forward or reverse).
        enzyme (Enzyme): Enzyme object containing the spacer length.
        direction (str, optional): The direction for spacer addition ("forward" or "reverse"). Default is "forward".

    Returns:
        str: A valid random spacer sequence.
    """
    spacer_length = enzyme.spacer_length

    while True:
        # Generate a random spacer of the required length
        spacer = generate_random_dna(spacer_length)

        # Add the spacer in the specified direction
        if direction == "forward":
            dna_new = rec_site + spacer + dna
        elif direction == "reverse":
            dna_new = dna + spacer + rec_site
        else:
            raise ValueError("Invalid direction. Use 'forward' or 'reverse'.")

        # Check if the new sequence passes the enzyme recognition site check
        if check_random_oligo_for_sites(dna_new, enzyme) == 1:
            return spacer



def generate_siteless_sequence(n, enzyme):
    """
    Generate a random DNA sequence that does not contain any enzyme recognition sites.

    Args:
        n (int): The length of the desired DNA sequence.
        enzyme (Enzyme): Enzyme object containing the recognition site properties.

    Returns:
        str: A random DNA sequence without any enzyme recognition sites.
    """
    while True:
        # Generate a random DNA sequence of the specified length
        random_oligo = generate_random_dna(n)

        # Check if the generated sequence contains any enzyme recognition sites
        if check_random_oligo_for_sites(random_oligo, enzyme) == 0:
            return random_oligo



def ensure_minimum_length(df, col, min_oligo_size, enzyme):
    """
    Ensure that each DNA sequence in the specified column of the DataFrame meets the minimum length requirement.

    Args:
        df (pd.DataFrame): Input DataFrame containing DNA sequences.
        col (str): The column name containing the DNA sequences to be processed.
        min_oligo_size (int): The minimum required length for each DNA sequence.
        enzyme (Enzyme): Enzyme object with recognition sequence properties.

    Returns:
        pd.DataFrame: DataFrame with updated DNA sequences that meet the minimum length.
    """
    for index, row in df.iterrows():
        modified_dna = row[col]

        # Ensure the modified DNA is a string (convert list/tuple to string if necessary)
        if not isinstance(modified_dna, str):
            if isinstance(modified_dna, (tuple, list)):
                modified_dna = ''.join(str(x) for x in modified_dna)
            else:
                modified_dna = str(modified_dna)

        current_length = len(modified_dna)

        # Calculate the number of bases to add (ensuring at least 6 bases are added)
        bases_to_add = max(np.floor((min_oligo_size - current_length) / 2), 6)

        # Generate random siteless sequences to meet length requirements
        random_bases1 = generate_siteless_sequence(bases_to_add, enzyme)
        random_bases2 = generate_siteless_sequence(bases_to_add, enzyme)

        # Concatenate random sequences with the original DNA sequence
        concatenated_DNA = random_bases1 + modified_dna + random_bases2

        df.at[index, col] = str(concatenated_DNA)

    return df
