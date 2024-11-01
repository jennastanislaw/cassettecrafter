import pandas as pd
from process_inputs_utils import (get_gene_name, get_dna_seq, 
                                get_allowed_codon_list, split_csl)

# Main function
def process_inputs(gene_file,mutations):
    name, starting_dna = read_input(gene_file)
    # TODO add some checks on the input in above fn

    mutations_df = mutation_file_to_df(mutations, starting_dna)

    return name, starting_dna, mutations_df


### Helper functions for process_inputs - sub-function read_input ###
def read_input(file):
    filetype=str(file).split("/")[-1].split(".")[-1]

    #check this is a file
    file_lines = open(file,"r").readlines()

    name = get_gene_name(file_lines,filetype)
    seq = get_dna_seq(file_lines,filetype)
    return name, seq

def mutation_file_to_df(mutations, og_seq_dna):
    og_seq_aa = og_seq_dna.translate()  # built-in biopython function
    mutation_df = pd.read_csv(mutations,index_col=0)

    # Add column that is list of mutations
    mutation_df["mut_list"] = mutation_df.iloc[:,0].apply(split_csl)

    og_aa = list()
    og_codon = list()
    for pos in mutation_df.index.tolist():
        pos_i = pos - 1 #adjsut pos because mutation indexing starts at 1, not 0
        og_aa.append(og_seq_aa[pos_i])
        og_codon.append("".join(og_seq_dna[3*pos_i:3*(pos_i+1)]))
    mutation_df["original"] = og_aa
    mutation_df["codons_original"] = og_codon

    mutation_df["allowed"] = mutation_df["mut_list"] + [mutation_df["original"]]

    # Add current codons and allowed codons (selecting the first one on the list)
    mutation_df["codons_allowed"] = get_allowed_codon_list(mutation_df["mut_list"].tolist(),
                                     mutation_df["codons_original"].tolist())

    return mutation_df