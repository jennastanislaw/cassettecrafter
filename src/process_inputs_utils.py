# import sys
import copy
from Bio.Seq import Seq
from dna_aa_definitions import CODON_TABLE_DNA, CODON_TO_AMINO_ACID_DNA, MIXED_BASES, MIXED_BASES_COMBO_TO_BASE

### Helper functions for read_input ###
# Get name from te fasta or csv
def get_gene_name(lines,filetype):
    """Get gene name from the first column of csv, or the first line of a fasta

    Args:
        lines (list): list of strings corresponding to lines in an input file
        filetype (str): string containing the file type of an input file 
            (determined by its extension)

    Returns:
        str : name as a string
    """
    # If csv, take name from the first column
    if filetype == "csv":
        name=lines[0].split(",")[0]
    # If fasta, take names from the first row
    elif filetype == "fa" or filetype == "fasta":
        name=lines[0].replace(">","")
    else:
        # return None
        raise ValueError(f"Filetype {filetype} is not allowed. Valid input file types are csv or fa")

    #remove extra characters if needed
    return name.strip().replace('\ufeff', '') 

# Read in the seq from the file
def get_seq(lines,filetype):
    """Get DNA sequence from the second column of csv, or the remaining lines (lines after name) of a fasta

    Args:
        lines (list): list of strings corresponding to lines in an input file
        filetype (str): string containing the file type of an input file 
            (determined by its extension)

    Returns:
        BioPython Seq: Seq object containing information about the original DNA 
            sequence of the gene
    """
    if filetype == "csv":
        seq=lines[0].split(",")[-1].strip()
    elif filetype == "fa":
        seq=""
        for line in lines[1:]:
            seq+=line.strip()

    #make it a biopython sequence for processing later
    # seq_obj = Seq(seq)

    # bases={"A","C","G","T"}

    # # If this is not a DNA sequence, convert it 
    # if not set(seq).issubset(bases):
    #     seq_obj = seq_obj.translate()
    
    return biopython_seq_from_str(seq)

def biopython_seq_from_str(str_seq):
    #make it a biopython sequence for processing later
    seq_obj = Seq(str_seq)

    bases={"A","C","G","T"}

    # If this is not a DNA sequence, convert it 
    if not set(str_seq).issubset(bases):
        seq_obj = convert_aa_to_dna(seq_obj._data)
    
    return seq_obj

def convert_aa_to_dna(prot_seq, codon_library=CODON_TABLE_DNA):
    dna_seq = ""
    allowed_aas = list(codon_library.keys())
    for aa in prot_seq.upper():
        if aa in allowed_aas:
            dna_seq += codon_library[aa][0] # select the first codon for the amino acid
        else:
            raise ValueError(f"Invalid amino acid: {aa}")
    
    return Seq(dna_seq)


### Helper functions for mutation_file_to_df ###
def gen_per_pos_muts(mutation_df):
    """Convert "allowed" column in dataframe to dictionary

    Args:
        mutation_df (pandas DataFrame): dataframe containing data about allowed
            mutations at different positions in the original sequence

    Returns:
        dict: dictionary mapping indices of input dataframe to values in the column 
            labelled "allowed"
    """
    # TODO: check that asserts df contains this column
    mutation_dict = mutation_df.loc["allowed"].to_dict()
    return mutation_dict

# Take comma separated list and return split list of strings
def split_csl(csl):
    """Splits an input string that contains commas into a list of string, 
        where each item in the list correponds a value in the original string list

        For example, for the input 'A,B,C', the output would be ['A','B','C']

    Args:
        csl (str): "comma-separated-list"; string containing a list that is delimited
            by commas

    Returns:
        list: list of strings
    """
    return list(val.strip() for val in csl.split(","))

def get_allowed_codon_list(mutations_aa, original_codons):
    """Convert amino acid information to codon information for each of the
        allowed mutations

    Args:
        mutations_aa (Pandas Series): columns from Pandas DataFrame which contains a list
            of allowed amino acid mutations for each row (position in the sequence)
        original_codons (Pandas Series): columns from Pandas DataFrame which contains 
            the original codon that the corresponding position

    Returns:
        list: list of allowed codons. The order of this list correlates with the
            order of positions given in the mutation csv and its dataframe. Note this
            includes the original codon sequence at that position, because it is
            also possible to have no mutation at a given position
    """
    output=list()

    # For each positions's possible mutations, convert the list of amino acids
    # to a list of codons
    for pos,allowed_aa_at_pos in enumerate(mutations_aa):
        allowed_codons_at_pos= list() #TODO: rename this variable
        original_codon = original_codons[pos]
        original_aa = CODON_TO_AMINO_ACID_DNA[original_codon]
    
        # Convert each AA in list to the first codon in the translation table
        for aa_mut in allowed_aa_at_pos:
            if not aa_mut == original_aa: #Don't add mutation to list if it is the same as the original
                allowed_codons_at_pos.append(CODON_TABLE_DNA[aa_mut][0])

        # Append list of allowed mutation codons and the original codon to output
        allowed_codons_at_pos.append(original_codon)

        final_allowed_codon_list = add_mixed_bases_and_combine(allowed_codons_at_pos)
        output.append(final_allowed_codon_list)

    return output

# fxn to combine any AA with shared first 2 bases
def add_mixed_bases_and_combine(allowed_codons_at_pos):
    modifiable_allowed_codon_list = copy.deepcopy(allowed_codons_at_pos)
    #TODO: add check of length of mutation list

    if len(allowed_codons_at_pos) < 1:
        return allowed_codons_at_pos
    for i, curr_codon in enumerate(allowed_codons_at_pos[:-1]):
        curr_codon_mod, base_1_2 = check_leu_arg_ser(curr_codon)
        # print("current", curr_codon)

        # Compare beginning of current codon to those in the allowed codons list
        for j, match_codon in enumerate(allowed_codons_at_pos[i+1:]):
            match_codon_mod, match_first_two = check_leu_arg_ser(match_codon)
            # print("match", match_codon)

            # If there is a match, find the mixed base that fits both codons
            if match_first_two == base_1_2: # add something for special cases: ARG (R) and LEU (L)
                mixed_base_codon = get_mixed_base_codon(curr_codon_mod, match_codon_mod)
                if mixed_base_codon[:2] != base_1_2:
                    mixed_base_codon = base_1_2 + mixed_base_codon[-1]

                # Remove the two codons and replace with the mixed-base codon
                modifiable_allowed_codon_list.remove(curr_codon)
                modifiable_allowed_codon_list.remove(match_codon)
                modifiable_allowed_codon_list.append(mixed_base_codon)

                break #exit the inner loop, because no other matches could exist
    return modifiable_allowed_codon_list 

def get_mixed_base_codon(codonA, codonB):
    base3_l = [codonA[-1],codonB[-1]]
    base3_l.sort()
    base3_joint_str = "".join(base3_l)

    # Access the 
    base3_mixed = MIXED_BASES_COMBO_TO_BASE[base3_joint_str]

    # Full codon with the third position as a mixed base
    base12 = codonA[:2]
    codon = base12 + base3_mixed

    return codon

def check_leu_arg_ser(codon):
    # Assign the first possible codon with a possible match to another amino acid codon if the first two bases
    codon_base_1_2 = codon[:2]
    if codon_base_1_2 == "CT": # Leu
        codon_base_1_2 = "TT"
        codon = "TTA"
    elif codon_base_1_2 == "CG": # Arg
        codon_base_1_2 = "AG"
        codon = "AGA"
    elif codon_base_1_2 == "TC": # Ser
        codon_base_1_2 = "AG"
        codon = "AGC"
    return codon, codon_base_1_2