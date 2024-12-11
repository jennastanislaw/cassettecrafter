"""
Processes and appropriately formats inputs to the program. The user must provide 
the sequence of the gene of interest, which can be given as either a path to a fasta
file or a csv file OR passed in as a string directly, as well as a csv
containing information about allowed/desired mutations that should be 
made to the sequence.

The gene sequence can be either an amino acid or a DNA sequence - the program
will convert a provided amino acid sequnce to DNA automatically, provided that
the protein sequence consists only of the 20 canonical amino acids.

An example of how the mutation csv should be formatted can be found at 
test_data/demo_mutation_list.csv. Briefly, this file should contain a column
called positions with the positions where mutations can occur, and a second
column called AminoAcids, which contains a list of allowed mutations. Thus, 
each row contain a position in the sequence an a list of allowed mutations at 
that postion.

After processing a tuple containing (the name of the sequence, the DNA sequence
as a Biopython Sequence object, the mutation information in a processed Pandas DataFrame)
"""

import pandas as pd
import os
from utils.process_inputs_utils import (read_gene_input_from_file,
                                        biopython_seq_from_str, 
                                        mutation_file_to_df)

def process_inputs(gene,mutations):
    """Process gene sequence input and allowed mutation input

    Args:
        gene (str): Either path to file containing genes to be inserted as fasta or csv
                    OR a string of the sequence. Can be a DNA or an amino acid sequence
        mutations (str): Path to csv file with mutation information

    Returns:
        tuple (name, starting_dna, mutations_df): tuple with gene name, initial 
                gene sequence, and pandas DataFrame with mutation information
    """
    # Check if gene input is a file or a string
    if os.path.isfile(gene):
        name, starting_dna_str = read_gene_input_from_file(gene)
    else:
        name = "gene1"
        starting_dna_str = gene
    starting_dna = biopython_seq_from_str(starting_dna_str)

    mutations_df = mutation_file_to_df(mutations, starting_dna)

    return name, starting_dna, mutations_df
