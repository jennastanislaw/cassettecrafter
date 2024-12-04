from Bio.Seq import Seq
import random
import numpy as np
from enzyme_site_replacement_utils import find_matching_sites, load_enzymes_from_csv, create_enzyme_dict


def generate_random_dna(n):

    """Generate a random DNA sequence of length n."""
    if n == 0:
        return Seq("")  # Return an empty Biopython Seq object for length 0

    if n < 0:
        raise ValueError("Length of DNA sequence must be a positive integer.")

    # Define the nucleotides
    nucleotides = ['A', 'T', 'C', 'G']

    # Generate a random sequence
    random_sequence = ''.join(random.choices(nucleotides, k=n))

    # Create a Biopython Seq object
    dna_sequence = Seq(random_sequence)

    return dna_sequence


def check_random_oligo_for_sites(DNA, enzyme):
    fwd_matches, rev_matches = find_matching_sites(enzyme, DNA)

    return len(rev_matches + fwd_matches)


# append end regions to all sequences
def add_end_restriction_sites(df, column, enzyme):
    """Add forward enzyme recognition site at the beginning and reverse recognition site at the end of each DNA sequence in the DataFrame."""

    # need to deal with spacer base length, also make sure that this variable is being used correctly everywhere else

    # Extract the forward and reverse recognition sites
    fwd_recognition_site = enzyme.fwd_recognition_site
    rev_recognition_site = enzyme.rev_recognition_site

    spacer_length = enzyme.spacer_length

    # Update each row in the DataFrame
    df[column] = df[column].apply(lambda dna: generate_spacers_rest_sites(dna, enzyme))

    return df


def generate_spacers_rest_sites(dna, enzyme):
    spacer_length = enzyme.spacer_length

    forward_recognition_site = enzyme.fwd_recognition_site

    reverse_recognition_site = enzyme.rev_recognition_site

    forward_rest_spacer = generate_spacer(dna, forward_recognition_site, enzyme, 'forward')

    reverse_rest_spacer = generate_spacer(dna, reverse_recognition_site, enzyme, 'reverse')

    dna_new = forward_recognition_site + forward_rest_spacer + dna + reverse_rest_spacer + reverse_recognition_site

    return str(dna_new)


def generate_spacer(dna, rec_site, enzyme, direction="forward"):
    """
    Generate a random spacer for a DNA sequence and ensure it passes the recognition site check.

    Parameters:
        dna (str): The DNA sequence.
        rec_site (str): The recognition site.
        enzyme (Enzyme): The enzyme object with recognition site properties.
        direction (str): Either "forward" or "reverse", indicating where to add the spacer.

    Returns:
        str: A valid random spacer sequence.
    """
    spacer_length = enzyme.spacer_length

    while True:

        spacer = generate_random_dna(spacer_length)

        # Add the spacer based on the direction
        if direction == "forward":
            dna_new = rec_site + spacer + dna
        elif direction == "reverse":
            dna_new = dna + spacer + rec_site
        else:
            raise ValueError("Invalid direction. Use 'forward' or 'reverse'.")

        if check_random_oligo_for_sites(dna_new, enzyme) == 1:
            return spacer


def generate_siteless_sequence(n, enzyme):
    while True:
        # Generate a random oligo of at least min_length
        random_oligo = generate_random_dna(n)  # Replace with your actual function to generate random oligos

        # Check if the generated random oligo has valid enzyme recognition sites
        if check_random_oligo_for_sites(random_oligo, enzyme) == 0:
            return random_oligo


def ensure_minimum_length(df, col, min_oligo_size, enzyme):
    """Ensure that each 'Modified DNA' sequence meets the minimum length requirement."""
    for index, row in df.iterrows():
        modified_dna = row[col]

        # Ensure `modified_dna` is a string
        if not isinstance(modified_dna, str):
            # Convert list or tuple elements to strings before joining
            if isinstance(modified_dna, (tuple, list)):
                modified_dna = ''.join(str(x) for x in modified_dna)
            else:
                modified_dna = str(modified_dna)

        current_length = len(modified_dna)

        # Calculate how many bases to add
        bases_to_add = max(np.floor((min_oligo_size - current_length) / 2), 6)  # Ensure at least 6 bases are added
        random_bases1 = generate_siteless_sequence(bases_to_add, enzyme)  # generate separate random sequences to meet
        # sequence complexity requirements

        random_bases2 = generate_siteless_sequence(bases_to_add, enzyme)
        # Append the random bases to the 'Modified DNA'

        concatenated_DNA = random_bases1 + modified_dna + random_bases2

        df.at[index, col] = str(concatenated_DNA)

    return df


