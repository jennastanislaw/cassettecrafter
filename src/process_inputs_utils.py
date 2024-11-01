from Bio.Seq import Seq
from dna_aa_definitions import CODON_TABLE_DNA

### Helper functions for read_input ###
# Get name from te fasta or csv
def get_gene_name(lines,filetype):
    if filetype == "csv":
        name=lines[0].split(",")[0]
    elif filetype == "fa":
        name=lines[0].replace(">","")
    return name.strip().replace('\ufeff', '') #remove extra characters if needed

# Read in the seq from the file
def get_dna_seq(lines,filetype):
    if filetype == "csv":
        seq=lines[0].split(",")[-1].strip()
    elif filetype == "fa":
        seq=""
        for line in lines[1:]:
            seq+=line.strip()

    return Seq(seq) #make it a biopython sequence for processing later


### Helper functions for mutation_file_to_df ###
def gen_per_pos_muts(mutation_df):
    mutation_dict = mutation_df.loc["allowed"].to_dict()
    return mutation_dict

# Take comma separated list and return split list of strings
def split_csl(csl):
    return [val.strip() for val in csl.split(",")]

def get_allowed_codon_list(mutations_aa, original_codons):
    output=list()
    for pos,allowed_aa_at_pos in enumerate(mutations_aa):
        output_pos= list()
        original_codon = original_codons[pos]
        for aa_mut in allowed_aa_at_pos:
            output_pos.append(CODON_TABLE_DNA[aa_mut][0])
        output_pos.append(original_codon)
        output.append(output_pos)
    return output