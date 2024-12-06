"""
From provided sequence and allowed mutations, produces a library of sequences
that contain all combinations of the desired mutations that are compatible with 
Golden Gate Assembly

Usage: python3 main.py -m [allowed mutations file] -b [plasmid backbone] 
        -e [restriction enzyme name] -d [file with enzyme data] 
        -m [minimum oligo size] -f [file with genes to insert]

Example: python3 src/main.py -m ./test_data/demo_mutation_list.csv 
                -f ./test_data/LY011_test_seq_single.csv
    When run from the directory above src. Adjust paths as needed

"""

import argparse
import os
import pandas as pd

from process_inputs import process_inputs
from generate_mutant_lib import generate_mutant_lib
from process_sequence_list import process_dna_sequences
from enzyme_site_replacement_utils import (load_enzymes_from_csv, create_enzyme_dict)
from mixed_base_rec_site_check import degen_codon_checker
from enzyme_site_replacement import remove_enzyme_sites
from split_sites import (find_split_indices)
from src.split_sites_utils import generate_cassettes


# to include in global variable file:
# enzyme_data
# Codon table to DNA
# mixed base codon info

# remove plasmid backbone

def generate_assembly_library(gene_file, mutations, enzyme_name, min_oligo_size, max_oligo_size):
    """Generates Golden Gate-compatible sequence library containing all possible
        combinations of allowed mutations 

    Args:
        gene_file (str): Path to file containing genes to be inserted. Can be a fasta or csv
        mutations (str): Path to csv file with mutation information
        enzyme_data (str): Path to file containing enzyme information
        enzyme_name (str): Name of enzyme whose recognition sites should be incorporated
                    into the genes. Note this name must match the name in the 
                    enzyme_data file
        min_oligo_size (int): Minimum oligo size required for each DNA sequence

    Returns:
        DataFrame : Pandas DataFrame containing mutation name and sequence
    """
    enzyme_fp = './data/enzyme_sites.csv'

# 1. Create enzyme object
    enzyme_class_objs = load_enzymes_from_csv(enzyme_fp) # should this be enzyme_data?
    enzyme_dict = create_enzyme_dict(enzyme_class_objs)

    enzyme = enzyme_dict.get(enzyme_name)

    if enzyme is None:
        raise KeyError(f"Enzyme '{enzyme_name}' not found.")

# 2. Process inputs
    name, starting_dna, mutations_df = process_inputs(gene_file,mutations)


# 3. Generate library of mutant sequences
    library_df = generate_mutant_lib(starting_dna,mutations_df, name)

    #this needs to be removed, along with pandas import
    library_df['DNA'] = (library_df['DNA'].str.replace(r"b'|b\"|'|\"", '', regex=True)  # Remove b', b", ', and "
    .str.replace(r'\s+', '', regex=True)          # Remove extra spaces, if any
)

# 4. Replace any unwanted enzyme sites
## This is going to make the naming of the sequences weird, will need to fix this

# 4.1 check degenerate codons for sites

    valid_mixed_bases = degen_codon_checker(library_df, enzyme)

# 4.2 replace sites

    rec_sites_removed = remove_enzyme_sites(enzyme, valid_mixed_bases)
# 5. Find split sites

    reference = rec_sites_removed['rec_sites_removed'].iloc[0]

    sequences = rec_sites_removed['rec_sites_removed'].tolist()

    split_sites = find_split_indices(reference, sequences, min_oligo_size, max_oligo_size, enzyme)

# 6 generate cassettes
    cassettes_df = generate_cassettes(rec_sites_removed, split_sites, enzyme)

    # 7. Final sequence processing (add terminal enzyme sites and extra bases if needed)

    final_df = process_dna_sequences(cassettes_df, enzyme, min_oligo_size)

    print(final_df) # export to csv in a way that pushes to the website

    return final_df

def parseargs():
    parser=argparse.ArgumentParser(
        description="""Generates Golden Gate-compatible sequence library containing all possible
        combinations of allowed mutations""", 
        prog="main.py"
    )

    parser.add_argument('--gene_file','-f', type=str,
                        default='',
                        help='File containing gene to insert')
    parser.add_argument('--mutations','-u', type=str,
                        default='', required=True,
                        help="""File containing mutations for genes. Should contain 
                        the amino acid position (starting from 1) in the first column,
                        and the allowed mutations in the second column""")
    parser.add_argument("--enzyme_name","-e",type=str,
                        default="BbsI", required=False,
                        help="""Enzyme name, matching name in the enzyme data file. Default is BbsI""")
    parser.add_argument("--min_oligo_size", "-m", type=int,
                        default=20, required=True,
                        help="""Minimum oligo size required for each DNA sequence. Default 20""")
    parser.add_argument("--max_oligo_size", "-M", type=int,
                        default=100, required=True,
                        help="""Maximum oligo size required for each DNA sequence. Default 100""")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args=parseargs()
    #Main function
    generate_assembly_library(args.gene_file, args.mutations, args.enzyme_name, args.min_oligo_size, args.max_oligo_size)


    # # Uncomment this line to run with hardcoded arguments for testing
    # final_result = generate_assembly_library(
    #     gene_file="/Users/siobhan/PycharmProjects/cassettecrafter/test_data/LY011_test_seq_single.csv",
    #     mutations="/Users/siobhan/PycharmProjects/cassettecrafter/test_data/demo_mutation_list.csv",
    #     enzyme_name="BbsI",
    #     min_oligo_size=20,
    #     max_oligo_size=100
    # )

    # # Print or process the final result if needed
    # print(final_result)