import pandas as pd
import os
from process_inputs_utils import (get_gene_name, get_dna_seq, 
                                get_allowed_codon_list, split_csl)

# Main function
def process_inputs(gene_file,mutations):
    """Process gene sequence input and allowed mutation input

    Args:
        gene_file (str): Path to file containing genes to be inserted. Can be a fasta or csv
        mutations (str): Path to csv file with mutation information

    Returns:
        tuple (name, starting_dna, mutations_df): tuple with gene name, initial 
                gene sequence, and pandas DataFrame with mutation information
    """
    name, starting_dna = read_dna_input(gene_file)
    # TODO add some checks on the input in above function

    mutations_df = mutation_file_to_df(mutations, starting_dna)

    return name, starting_dna, mutations_df


### Helper functions for process_inputs - sub-function read_input ###
def read_dna_input(file):
    """Extract gene name and sequence from input file

    Args:
        file (str):  Path to file containing genes to be inserted. Can be a fasta or csv

    Returns:
        tuple (name, seq): tuple with gene name and initial gene sequence (DNA)
    """
    filetype=str(file).split("/")[-1].split(".")[-1]

    #check this is a file
    file_lines = open(file,"r").readlines()

    # TODO: add functionality to accept amino acide sequence instead of DNA
    name = get_gene_name(file_lines,filetype)
    seq = get_dna_seq(file_lines,filetype)
    return name, seq

def mutation_file_to_df(mutations, og_seq_dna):
    """Generates a dataframe with information about the allowed mutations

    Args:
        mutations (str): Path to csv file with mutation information
        og_seq_dna (BioPython Seq): BioPython Seq object containing information
            about the original DNA sequence of the gene

    Returns:
        pandas DataFrame : dataframe containing data from the original csv and
            additional columns containng details about each mutatable position: 
            [original - original amino acide, codons_original - original codon,
            allowed - list of allowed amino acids at that position,
            codons_allowed - list of allowed codons at that position (corresponds 
            to allowed amino acids)]
    """
    og_seq_aa = og_seq_dna.translate()  # built-in biopython function
    mutation_df = pd.read_csv(mutations,index_col=0)

    # Add column that is list of mutations
    # TODO: modify this so mutation information can be formatted differently
    #   (e.g. without commas, or as a space-separated list)
    mutation_df["mut_list"] = mutation_df.iloc[:,0].apply(split_csl)

    og_aa = list()
    og_codon = list()
    for pos in mutation_df.index.tolist():
        pos_i = pos - 1 #adjsut pos because mutation indexing starts at 1, not 0
        og_aa.append(og_seq_aa[pos_i])
        og_codon.append("".join(og_seq_dna[3*pos_i:3*(pos_i+1)]))
    mutation_df["original"] = og_aa
    mutation_df["codons_original"] = og_codon

    # mutation_df["allowed"] = mutation_df["mut_list"] + [mutation_df["original"]]
    mutation_df["allowed"] = mutation_df.apply(lambda row: row['mut_list'] + [row['original']], axis=1)

    # Add current codons and allowed codons (selecting the first one on the list)
    mutation_df["codons_allowed"] = get_allowed_codon_list(mutation_df["mut_list"].tolist(),
                                     mutation_df["codons_original"].tolist())
    
    print( mutation_df)

    return mutation_df