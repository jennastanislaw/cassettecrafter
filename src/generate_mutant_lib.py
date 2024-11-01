import pandas as pd
import itertools
from generate_mutant_lib_utils import split_to_codons, make_mut_dict


def generate_mutant_lib(og_codon_seq,mutations_df, name):
    # Load mutation file and generate dictionary of allowed mutations
    mutation_dict = mutations_df["codons_allowed"].to_dict()

    # Generate all combinations of possible mutations (including the original codons)
    all_combinations = list(itertools.product(*[value for key, value in 
                                                mutation_dict.items()]))
    # Split dna sequence into codons - group every 3 bases
    editable_codon_seq=split_to_codons(og_codon_seq._data)

    # Generate dictionary where keys are names for mutated sequence and 
    # values are the codon sequence corresponding to that mutant
    mut_library = make_mut_dict(editable_codon_seq, all_combinations, name, mutations_df)

    # Turn dict into dataframe for next steps
    indices = list(mut_library.keys())
    seq_col = list(mut_library.values())
    library_df = pd.DataFrame({"name":indices,"DNA": seq_col}) 

    return library_df
