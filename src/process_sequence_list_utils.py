from Bio.Seq import Seq
import random
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