import pandas as pd
from Bio.Seq import Seq
import itertools
from enzyme_site_replacement_utils import codon_table_dna,get_codons_with_full_indices # should capitalize this -global 

def read_input(file):
    filetype=str(file).split("/")[-1].split(".")[-1]

    #check this is a file
    file_lines = open(file,"r").readlines()

    name = get_gene_name(file_lines,filetype)
    seq = get_dna_seq(file_lines,filetype)
    return name, seq

def get_gene_name(lines,filetype):
    if filetype == "csv":
        name=lines[0].split(",")[0]
    elif filetype == "fa":
        name=lines[0].replace(">","")
    return name.strip().replace('\ufeff', '')

# Read in the sequnence from the file
def get_dna_seq(lines,filetype):
    if filetype == "csv":
        seq=lines[0].split(",")[-1].strip()
    elif filetype == "fa":
        seq=""
        for line in lines[1:]:
            seq+=line.strip()

    return Seq(seq) #make it a biopython sequence for processing later

# Take comma separated list and return split list of strings
def split(csl):
    return [val.strip() for val in csl.split(",")]

def get_allowed_codon_list(mutations_aa, original_codons):
    output=list()
    for pos,allowed_aa_at_pos in enumerate(mutations_aa):
        output_pos= list()
        original_codon = original_codons[pos]
        for aa_mut in allowed_aa_at_pos:
            output_pos.append(codon_table_dna[aa_mut][0])
        output_pos.append(original_codon)
        output.append(output_pos)
    return output


def mutation_file_to_df(mutations, og_seq_dna):
    og_seq_aa = og_seq_dna.translate()  # built-in biopython function
    mutation_df = pd.read_csv(mutations,index_col=0)

    # Add column that is list of mutations
    mutation_df["mut_list"] = mutation_df.iloc[:,0].apply(split)

    og_aa = list()
    og_codon = list()
    for pos in mutation_df.index.tolist():
        pos_i = pos - 1 #adjsut pos because mutation indexing starts at 1, not 0
        og_aa.append([og_seq_aa[pos_i]])
        og_codon.append("".join(og_seq_dna[3*pos_i:3*(pos_i+1)]))
    mutation_df["original"] = og_aa
    mutation_df["codons_original"] = og_codon

    mutation_df["allowed"] = mutation_df["mut_list"] + mutation_df["original"]

    # Add current codons and allowed codons (selecting the first one on the list)
    mutation_df["codons_allowed"] = get_allowed_codon_list(mutation_df["mut_list"].tolist(),
                                     mutation_df["codons_original"].tolist())

    return mutation_df


def gen_per_pos_muts(mutation_df):
    mutation_dict = mutation_df.loc["allowed"].to_dict()
    return mutation_dict

def generate_mutant_lib(og_codon_seq,mutations_df, name):
    # Load mutation file and generate dictionary of allowed mutations
    # mutation_df = mutation_file_to_df(mutations, og_seq)
    mutation_dict = mutations_df["allowed"].to_dict()

    # Generate all combinations of possible mutations (including the original codons)
    all_combinations = list(itertools.product(*[value for key, value in 
                                                mutation_dict.items()]))

    positions=mutations_df.index.tolist()
    editable_codon_seq=get_codons_with_full_indices(og_codon_seq)
    # editable_seq = list(og_seq)
    
    mut_library = dict()
    for combo in all_combinations:
        mut_name=name
        for i,pos in enumerate(positions):
            editable_codon_seq[pos-1] = combo[i]
            mut_name+=f"_{str(pos)}{combo[i]}"
        seq="".join(editable_codon_seq)
        mut_library[mut_name] = seq

    return mut_library
