from Bio.Seq import Seq
import random
import numpy as np
from enzyme_site_replacement_utils import find_matching_sites, load_enzymes_from_csv, create_enzyme_dict
def generate_random_dna(n):
    """Generate a random DNA sequence of length n."""
    if n <= 0:
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

    if len(rev_matches+fwd_matches) > 0 :
        ok = False
    else:
        ok = True
    return ok


# append end regions to all sequences
def add_end_restriction_sites(df, enzyme_fp, enzyme, min_oligo_size):

    """Add forward enzyme recognition site at the beginning and reverse recognition site at the end of each DNA sequence in the DataFrame."""

    # Extract the forward and reverse recognition sites
    fwd_recognition_site = enzyme.fwd_recognition_site
    rev_recognition_site = enzyme.rev_recognition_site

    # Ensure that the 'DNA' column exists in the DataFrame
    if 'DNA' not in df.columns:
        raise ValueError("DataFrame must contain a 'DNA' column.")

    # Update each row in the DataFrame
    df['Modified_DNA'] = df['DNA'].apply(lambda dna: fwd_recognition_site + dna + rev_recognition_site)

    return df

def generate_siteless_sequence(n, enzyme):
    while True:
        # Generate a random oligo of at least min_length
        random_oligo = generate_random_dna(n)  # Replace with your actual function to generate random oligos

        # Check if the generated random oligo has valid enzyme recognition sites
        if check_random_oligo_for_sites(random_oligo, enzyme):
            return random_oligo

def ensure_minimum_length(df, min_oligo_size, enzyme_name):
    """Ensure that each 'Modified DNA' sequence meets the minimum length requirement."""
    for index, row in df.iterrows():
        modified_dna = row['Modified_DNA']
        current_length = len(modified_dna)

        # Calculate how many bases to add
        bases_to_add = max(np.floor((min_oligo_size - current_length)/2), 6)  # Ensure at least 6 bases are added
        random_bases1 = generate_siteless_sequence(bases_to_add, enzyme_name) # generate separate random sequences to meet
        # sequence complexity requirements

        random_bases2 = generate_siteless_sequence(bases_to_add, enzyme_name)
        # Append the random bases to the 'Modified DNA'
        df.at[index, 'Modified_DNA'] = random_bases1 + modified_dna + random_bases2

    return df

