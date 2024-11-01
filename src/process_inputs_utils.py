import sys
from Bio.Seq import Seq
from dna_aa_definitions import CODON_TABLE_DNA

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
        raise ValueError(f"Filetype {filetype} is not allowed. Valid input file types are csv or fa")

    #remove extra characters if needed
    return name.strip().replace('\ufeff', '') 

# Read in the seq from the file
def get_dna_seq(lines,filetype):
    """Get gene name from the first column of csv, or the first line of a fasta

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

    return Seq(seq) #make it a biopython sequence for processing later


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
    return [val.strip() for val in csl.split(",")]

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
    # to a list of codon
    for pos,allowed_aa_at_pos in enumerate(mutations_aa):
        output_pos= list()
    
        # Convert each AA in list to the first codon in the translation table
        for aa_mut in allowed_aa_at_pos:
            output_pos.append(CODON_TABLE_DNA[aa_mut][0])
        
        # Append list of allowed mutation codons and the original codon to output
        original_codon = original_codons[pos]
        output_pos.append(original_codon)
        output.append(output_pos)

    return output