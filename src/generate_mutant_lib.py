"""
From dataframe containing information about allowed mutations and the original
amino acids at those positions, generate codon sequences of all possible combinations
of mutations.

For example, in a gene with amino acid sequence MTLIVLHH, with allowed mutations
to position 2 of H and W, and allowed mutations to position 4 of L, the following
combination of amino acid sequences are allowed:
MTLIVLHH (original)
MHLIVLHH (T2H)
MHLLVLHH (T2H_I5L)
MWLIVLHH (T2W)
MWLLVLHH (T2W_I5L)
MTLLVLHH (I4L)

In this case, each of the above possible sequences would be converted to a DNA
sequence, and the output of the function will be a dataframe with the name of 
the mutant (as shown in parentheses above) linked to the DNA sequence of that mutant

"""

import pandas as pd
import itertools
from utils.generate_mutant_lib_utils import split_to_codons, make_mut_dict


def generate_mutant_lib(og_codon_seq,mutations_df, name):
    """From dataframe containing information about allowed mutations and the original
        amino acids at those positions, generate codon sequences of all possible combinations
        of mutations.

    Args:
        og_codon_seq (BioPython Seq): Seq object containing information about 
            the original DNA sequence of the gene
        mutations_df (pandas DataFrame): dataframe containing data about allowed
            mutations at different positions in the original sequence
        name (str): name of the gene whose sequences are being generates

    Returns:
        pandas DataFrame : dataframe containing mutant names and a "DNA" column
            with the DNA sequence for that mutant
    """
    # Load mutation file and generate dictionary of allowed mutations
    mutation_dict = mutations_df["codons_allowed"].to_dict()

    # Generate all combinations of possible mutations (including the original codons)
    all_combinations = list(itertools.product(*[value for key, value in 
                                                mutation_dict.items()]))

    # Split dna sequence into codons - group every 3 bases
    #._data to get sequence from BioPython Seq object
    editable_codon_seq=split_to_codons(og_codon_seq._data) 

    # Generate dictionary where keys are names for mutated sequence and 
    # values are the codon sequence corresponding to that mutant
    mut_library = make_mut_dict(editable_codon_seq, all_combinations, name, mutations_df)

    # Turn dict into dataframe for next steps
    indices = list(mut_library.keys())
    seq_col = list(mut_library.values())
    library_df = pd.DataFrame({"name":indices,"DNA": seq_col}) 

    return library_df
