import pandas as pd
from Bio.Seq import Seq
from dna_aa_definitions import CODON_TABLE_DNA, CODON_TO_AMINO_ACID_DNA

class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = int(spacer_length)
        self.OH_length = int(OH_length)

    def __repr__(self):
        return (f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, "
                f"rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, "
                f"OH_length={self.OH_length})")


# Function to create a dictionary of enzymes
def create_enzyme_dict(enzymes):
    """Creates a dictionary with enzyme names as keys and Enzyme objects as values."""
    return {enzyme.name: enzyme for enzyme in enzymes}

def load_enzymes_from_csv(csv_file_path):
    """Reads a CSV file and creates Enzyme objects for each row."""
    # Read the CSV file using pandas
    df = pd.read_csv(csv_file_path)

    enzymes = []
    # Iterate through each row in the DataFrame and create Enzyme objects
    for _, row in df.iterrows():
        enzyme = Enzyme(
            name=row['Enzyme'],
            fwd_recognition_site=row['Fwd_recognition_site'],
            rev_recognition_site=row['Rev_recognition_site'],
            spacer_length=row['spacer_length'],
            OH_length=row['OH_length']
        )
        enzymes.append(enzyme)
    return enzymes


#def find_IIs_sites(DNA, enzyme):
def find_matching_sites(enzyme, dna_sequence):
    """Finds positions of the forward and reverse recognition sites in a DNA sequence."""
    fwd_matches = []
    rev_matches = []

    # Find all positions of the forward recognition site
    fwd_site_length = len(enzyme.fwd_recognition_site)
    rev_site_length = len(enzyme.rev_recognition_site)

    for i in range(len(dna_sequence) - fwd_site_length + 1):
        if dna_sequence[i:i + fwd_site_length] == enzyme.fwd_recognition_site:
            fwd_matches.append(i)

    # Find all positions of the reverse recognition site
    for i in range(len(dna_sequence) - rev_site_length + 1):
        if dna_sequence[i:i + rev_site_length] == enzyme.rev_recognition_site:
            rev_matches.append(i)

    return fwd_matches, rev_matches

def generate_synonymous_codons_dna(codon):
    """Generate synonymous codons (alternative codons for the same amino acid) for DNA."""
    if codon not in CODON_TO_AMINO_ACID_DNA:
        raise ValueError(f"Invalid codon: {codon}")

    amino_acid = CODON_TO_AMINO_ACID_DNA[codon]
    synonymous_codons = CODON_TABLE_DNA[amino_acid]

    # Remove the input codon from the list to get alternatives
    return [c for c in synonymous_codons if c != codon]


def get_codons_with_full_indices(dna_sequence):
    """Returns a list of codons with their base indices in the DNA sequence."""
    seq = Seq(dna_sequence)
    codon_list = []

    # Iterate over the sequence in steps of 3 to get codons
    for i in range(0, len(seq), 3):
        codon = str(seq[i:i + 3])  # Extract the codon
        # Create a list of indices for the bases within this codon
        indices = [i, i + 1, i + 2]  # Indices for the three bases
        codon_list.append((codon, indices))  # Append a tuple of (codon, indices)

    return codon_list


def get_affected_codons_by_recognition_sites(dna_sequence, recognition_idxs, recognition_site):
    """Returns a dictionary with recognition site indices as keys and lists of affected codons and their indices as values."""
    codons_with_indices = get_codons_with_full_indices(dna_sequence)
    affected_codons_dict = {}

    for recognition_idx in recognition_idxs:
        recognition_site_indices = list(range(recognition_idx, recognition_idx + len(recognition_site)))

        # Find affected codons for the current recognition site
        affected_codons = []
        for codon, indices in codons_with_indices:
            if any(index in indices for index in recognition_site_indices):
                affected_codons.append((codon, indices))  # Store codon and its indices

        # Store in the dictionary
        affected_codons_dict[recognition_idx] = affected_codons

    return affected_codons_dict

