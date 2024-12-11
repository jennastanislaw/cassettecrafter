"""
This script provides functionality for working with enzymes, recognition sites, and DNA sequences. It includes a class
for representing enzymes with their recognition sites and properties, as well as several functions for analyzing DNA
sequences and identifying codons affected by recognition sites.

Classes:
    Enzyme: A class representing an enzyme with its recognition sites and properties, including the enzyme name, forward
    and reverse recognition sites, spacer length, and overhang length.

Functions:
    create_enzyme_dict: Creates a dictionary mapping enzyme names to their corresponding Enzyme objects.
    load_enzymes_from_csv: Reads a CSV file and creates a list of Enzyme objects from each row.
    find_matching_sites: Finds the positions of the forward and reverse recognition sites for an enzyme in a DNA
    sequence.
    generate_synonymous_codons_dna: Generates synonymous codons (alternative codons for the same amino acid) for a given
    DNA codon.
    get_codons_with_full_indices: Returns a list of codons along with their base indices in the DNA sequence.
    get_affected_codons_by_recognition_sites: Returns a dictionary with recognition site indices as keys and lists of
    affected codons and their indices as values.

Dependencies:
    pandas: For handling CSV files and creating DataFrames.
    Bio.Seq: For DNA sequence manipulation and codon extraction.
    dna_aa_definitions: For codon-to-amino acid mapping (CODON_TABLE_DNA, CODON_TO_AMINO_ACID_DNA).
"""
import pandas as pd
from Bio.Seq import Seq
from data.dna_aa_definitions import CODON_TABLE_DNA, CODON_TO_AMINO_ACID_DNA


class Enzyme:
    """
    A class representing an enzyme with its recognition sites and properties.

    Attributes:
        name (str): The name of the enzyme.
        fwd_recognition_site (str): The forward recognition site sequence for the enzyme.
        rev_recognition_site (str): The reverse recognition site sequence for the enzyme.
        spacer_length (int): The spacer length between recognition site and cut site.
        OH_length (int): The length of overhangs created by the enzyme.
    """

    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        """
        Initializes an Enzyme object with the given properties.

        Args:
            name (str): The name of the enzyme.
            fwd_recognition_site (str): The forward recognition site sequence.
            rev_recognition_site (str): The reverse recognition site sequence.
            spacer_length (int): The spacer length between recognition site and cut site.
            OH_length (int): The length of overhangs created by the enzyme.
        """
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = int(spacer_length)
        self.OH_length = int(OH_length)

    def __repr__(self):
        """
        Provides a string representation of the Enzyme object.

        Returns:
            str: A string describing the enzyme and its properties.
        """
        return (f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, "
                f"rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, "
                f"OH_length={self.OH_length})")


def create_enzyme_dict(enzymes):
    """
    Creates a dictionary mapping enzyme names to their corresponding Enzyme objects.

    Args:
        enzymes (list): A list of Enzyme objects.

    Returns:
        dict: A dictionary where keys are enzyme names and values are Enzyme objects.
    """
    # Use a dictionary comprehension to map enzyme names to Enzyme objects.
    return {enzyme.name: enzyme for enzyme in enzymes}


def load_enzymes_from_csv(csv_file_path):
    """
    Reads a CSV file and creates a list of Enzyme objects for each row.

    Args:
        csv_file_path (str): The path to the CSV file containing enzyme data.

    Returns:
        list: A list of Enzyme objects created from the CSV file.
    """
    # Read the CSV file into a pandas DataFrame.
    df = pd.read_csv(csv_file_path)

    enzymes = []
    # Iterate through each row in the DataFrame to create Enzyme objects.
    for _, row in df.iterrows():
        enzyme = Enzyme(
            name=row['Enzyme'],  # Enzyme name.
            fwd_recognition_site=row['Fwd_recognition_site'],  # Forward recognition site sequence.
            rev_recognition_site=row['Rev_recognition_site'],  # Reverse recognition site sequence.
            spacer_length=row['spacer_length'],  # Spacer length between recognition sites.
            OH_length=row['OH_length']  # Overhang length created by the enzyme.
        )
        enzymes.append(enzyme)

    return enzymes


def find_matching_sites(enzyme, dna_sequence):
    """
    Finds the positions of the forward and reverse recognition sites for an enzyme in a DNA sequence.

    Args:
        enzyme (Enzyme): An Enzyme object containing forward and reverse recognition sites.
        dna_sequence (str): The DNA sequence to search.

    Returns:
        tuple: Two lists containing the start positions of the forward and reverse recognition sites, respectively.
    """
    fwd_matches = []  # List to store forward recognition site positions.
    rev_matches = []  # List to store reverse recognition site positions.

    # Calculate the lengths of the recognition sites.
    fwd_site_length = len(enzyme.fwd_recognition_site)
    rev_site_length = len(enzyme.rev_recognition_site)

    # Find all positions of the forward recognition site in the DNA sequence.
    for i in range(len(dna_sequence) - fwd_site_length + 1):
        if dna_sequence[i:i + fwd_site_length] == enzyme.fwd_recognition_site:
            fwd_matches.append(i)

    # Find all positions of the reverse recognition site in the DNA sequence.
    for i in range(len(dna_sequence) - rev_site_length + 1):
        if dna_sequence[i:i + rev_site_length] == enzyme.rev_recognition_site:
            rev_matches.append(i)

    # Return the lists of forward and reverse site positions.
    return fwd_matches, rev_matches


def generate_synonymous_codons_dna(codon):
    """
    Generates synonymous codons (alternative codons for the same amino acid) for a given DNA codon.

    Args:
        codon (str): The DNA codon for which synonymous codons are to be generated.

    Returns:
        list: A list of synonymous codons for the same amino acid, excluding the input codon.

    Raises:
        ValueError: If the input codon is not valid (not found in the codon-to-amino acid mapping).
    """
    # Check if the codon is valid by looking it up in the codon-to-amino acid mapping.
    if codon not in CODON_TO_AMINO_ACID_DNA:
        raise ValueError(f"Invalid codon: {codon}")

    # Retrieve the amino acid corresponding to the given codon.
    amino_acid = CODON_TO_AMINO_ACID_DNA[codon]

    # Get the list of synonymous codons for the amino acid.
    synonymous_codons = CODON_TABLE_DNA[amino_acid]

    # Return the list of synonymous codons, excluding the original codon.
    return [c for c in synonymous_codons if c != codon]


def get_codons_with_full_indices(dna_sequence):
    """
    Returns a list of codons along with their base indices in the DNA sequence.

    Args:
        dna_sequence (str): The DNA sequence to extract codons from.

    Returns:
        list: A list of tuples, where each tuple contains a codon and a list of its base indices.

    Example:
        get_codons_with_full_indices("ATGCGT")
        returns [('ATG', [0, 1, 2]), ('CGT', [3, 4, 5])]
    """
    seq = Seq(dna_sequence)  # Convert the DNA sequence to a Seq object for easier slicing
    codon_list = []  # Initialize an empty list to store codons with indices

    # Iterate over the sequence in steps of 3 to extract codons
    for i in range(0, len(seq), 3):
        codon = str(seq[i:i + 3])  # Extract a 3-base codon
        indices = [i, i + 1, i + 2]  # List of indices for the three bases in the codon
        codon_list.append((codon, indices))  # Append the codon and its indices as a tuple

    return codon_list


def get_affected_codons_by_recognition_sites(dna_sequence, recognition_idxs, recognition_site):
    """
    Returns a dictionary with recognition site indices as keys and lists of affected codons and their indices as values.

    Args:
        dna_sequence (str): The DNA sequence to analyze.
        recognition_idxs (list): A list of indices where the recognition site starts.
        recognition_site (str): The recognition site sequence.

    Returns:
        dict: A dictionary where keys are recognition site indices, and values are lists of affected codons and their
        base indices.

    Example:
        get_affected_codons_by_recognition_sites("ATGCGT", [0], "ATGC")
        returns {0: [('ATG', [0, 1, 2]), ('CGT', [3, 4, 5])]}
    """
    codons_with_indices = get_codons_with_full_indices(dna_sequence)  # Get all codons with their indices
    affected_codons_dict = {}  # Initialize an empty dictionary to store affected codons

    # Iterate through each recognition site index
    for recognition_idx in recognition_idxs:
        recognition_site_indices = list(
            range(recognition_idx, recognition_idx + len(recognition_site)))  # Indices of the recognition site

        # Find all codons affected by the current recognition site
        affected_codons = []
        for codon, indices in codons_with_indices:
            if any(index in indices for index in
                   recognition_site_indices):  # Check if any base index is part of the recognition site
                affected_codons.append((codon, indices))  # Add the affected codon and its indices

        # Store the affected codons in the dictionary with the recognition site index as the key
        affected_codons_dict[recognition_idx] = affected_codons

    return affected_codons_dict

