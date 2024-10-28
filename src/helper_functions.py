import pandas as pd
from Bio.Seq import Seq
import itertools

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

def mutation_file_to_dict(mutations, og_seq):
    mutation_df = pd.read_csv(mutations,index_col=0)
    mutation_df["mut_list"] = mutation_df.iloc[:,0].apply(split)

    og = list()
    for pos in mutation_df.index.tolist():
        og.append([og_seq[pos-1]])
    mutation_df["original"] = og

    # print(mutation_df["mut_list"], mutation_df["original"])

    mutation_df["allowed"] = mutation_df["mut_list"] + mutation_df["original"]

    return mutation_df
    # mutation_dict = mutation_df["allowed"].to_dict()

    # return mutation_dict
    

def gen_per_pos_muts(mutation_df):
    mutation_dict = mutation_df.loc["allowed"].to_dict()
    return mutation_dict

def generate_mutant_lib(og_seq,mutations, name):
    mutation_df = mutation_file_to_dict(mutations, og_seq)
    mutation_dict = mutation_df["allowed"].to_dict()

    all_combinations = list(itertools.product(*[value for key, value in 
                                                mutation_dict.items()]))

    positions=mutation_df.index.tolist()
    editable_seq = list(og_seq)
    
    mut_library = dict()
    for combo in all_combinations:
        mut_name=name
        for i,pos in enumerate(positions):
            editable_seq[pos-1] = combo[i]
            mut_name+=f"_{str(pos)}{combo[i]}"
        seq="".join(editable_seq)
        mut_library[mut_name] = seq

    return mut_library
