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

from process_inputs import process_inputs
from generate_mutant_lib import generate_mutant_lib
from process_sequence_list import (replace_enzyme_sites_in_dataframe, 
                                    process_dna_sequences)
from enzyme_site_replacement_utils import (load_enzymes_from_csv, create_enzyme_dict)
from mixed_base_rec_site_check import degen_codon_checker
from enzyme_site_replacement import remove_enzyme_sites
from split_sites import (find_split_indices, generate_cassettes)

# to include in global variable file:
# enzyme_data
# Codon table to DNA
# mixed base codon info

# remove plasmid backbone

def generate_assembly_library(gene_file, mutations, backbone, enzyme_data, 
                      enzyme_name, min_oligo_size, max_oligo_size): # need to add in user specified min and max oligo lengths
    """Generates Golden Gate-compatible sequence library containing all possible
        combinations of allowed mutations 

    Args:
        gene_file (str): Path to file containing genes to be inserted. Can be a fasta or csv
        mutations (str): Path to csv file with mutation information
        backbone (str): Path to file containing sequence of the plasmid backbone
                    into which genes will be inserted. NOTE: this is currently not 
                    being utilized, but will be functional in future versions
        enzyme_data (str): Path to file containing enzyme information
        enzyme_name (str): Name of enzyme whose recognition sites should be incorporated
                    into the genes. Note this name must match the name in the 
                    enzyme_data file
        min_oligo_size (int): Minimum oligo size required for each DNA sequence

    Returns:
        DataFrame : Pandas DataFrame containing mutation name and sequence
    """

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
    cassettes_df = generate_cassettes(rec_sites_removed, split_sites)

    # 7. Final sequence processing (add terminal enzyme sites and extra bases if needed)

    final_df = process_dna_sequences(cassettes_df, enzyme, min_oligo_size)

    # print(final_df)

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
    parser.add_argument('--backbone','-b', type=str,
                        default='',required=False,
                        help='File containing plasmid backbone sequence. Note should contain RE sites for GG')
    parser.add_argument('--mutations','-m', type=str,
                        default='', required=False,
                        help="""File containing mutations for genes. Should contain 
                        the amino acid position (starting from 1) in the first column,
                        and the allowed mutations in the second column""")
    parser.add_argument("--enzyme_data", "-d", type=str, required=False,
                        default=f"{os.path.dirname(os.path.abspath(__file__))}/data/enzyme_sites.csv", 
                        help="""File path to enzyme data file. If no file is provided, the
                            default path will be used, referencing the provided enzyme data file
                            at src/data/enzyme_sites.csv""")
    parser.add_argument("--enzyme_name","-e",type=str,
                        default="BbsI", required=False,
                        help="""Enzyme name, matching name in the enzyme data file. Default is BbsI""")
    parser.add_argument("--min_oligo_size", "-s", type=int,
                        default=20, required=False,
                        help="""Minimum oligo size required for each DNA sequence. Default 20""")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args=parseargs()
    # Main function
    generate_assembly_library(args.gene_file, args.mutations, args.backbone,
                      args.enzyme_data, args.enzyme_name, args.min_oligo_size) 