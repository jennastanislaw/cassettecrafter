"""
From provided sequence and allowed mutations, produces a library of sequences
that contain all combinations of the desired mutations that are compatible with 
Golden Gate Assembly

Usage: python3 main.py -f [gene to insert] -u [allowed mutations file]
        -e [restriction enzyme name] -m [minimum oligo size]  -M [minimum oligo size]

Example: python3 src/main.py -f ./test_data/LY011_test_seq_single.csv
    -u ./test_data/demo_mutation_list.csv -M 50 -m 30

"""

import argparse
import os
import pandas as pd

from process_inputs import process_inputs
from generate_mutant_lib import generate_mutant_lib
from process_sequence_list import process_dna_sequences
from utils.enzyme_site_replacement_utils import (load_enzymes_from_csv, create_enzyme_dict)
from mixed_base_rec_site_check import degen_codon_checker
from enzyme_site_replacement import remove_enzyme_sites
from split_sites import (find_split_indices)
from utils.split_sites_utils import generate_cassettes
from process_output import process_output


def generate_assembly_library(gene, mutations, enzyme_name, min_oligo_size, max_oligo_size, output=""):
    """Generates Golden Gate-compatible sequence library containing all possible
        combinations of allowed mutations 

    Args:
        gene (str): Either path to file containing genes to be inserted as fasta or csv
                    OR a string of the sequence. Can be a DNA or an amino acid sequence
        mutations (str): Path to csv file with mutation information
        enzyme_name (str): Name of enzyme whose recognition sites should be incorporated
                    into the genes. Note this name must match the name in the 
                    enzyme_data file (located at src/data/enzyme_sites.csv)
        min_oligo_size (int): Minimum oligo size required for each DNA sequence
        max_oligo_size (int): Maximum oligo size required for each DNA sequence
        output (str): Path where output csv should be dumped.

    Returns:
        DataFrame : Pandas DataFrame containing mutation name and sequence
    """
    enzyme_data = f'{os.path.dirname(__file__)}/data/enzyme_sites.csv'
    # 1. Create enzyme object
    # enzyme_class_objs = load_enzymes_from_csv(enzyme_fp) # should this be enzyme_data? - yes
    enzyme_class_objs = load_enzymes_from_csv(enzyme_data)
    enzyme_dict = create_enzyme_dict(enzyme_class_objs)
    print(enzyme_dict)
    enzyme = enzyme_dict.get(enzyme_name)

    if enzyme is None:
        raise KeyError(f"Enzyme '{enzyme_name}' not found.")

    # 2. Process inputs
    name, starting_dna, mutations_df = process_inputs(gene,mutations)

    # 3. Generate library of mutant sequences
    library_df = generate_mutant_lib(starting_dna,mutations_df, name)

    # 4. Replace any unwanted enzyme sites

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
    filtered_df = process_output(final_df, output)

    return filtered_df

def parseargs():
    parser=argparse.ArgumentParser(
        description="""Generates Golden Gate-compatible sequence library containing all possible
        combinations of allowed mutations""", 
        prog="main.py"
    )

    parser.add_argument('--gene_file','-f', type=str,
                        default='',
                        help="""Sequence of gene to mutate for Golden Gate. Usually 
                                provided as a fasta file or a csv but can also be a string.
                                Note that the sequence can be DNA or protein, but 
                                DNA is recommended as this program does not allow
                                for organism-specific codon optimization""")
    parser.add_argument('--mutations','-u', type=str,
                        default='', required=True,
                        help="""File containing mutations for genes. Should contain 
                        the amino acid position (starting from 1) in the first column,
                        and the allowed mutations in the second column""")
    parser.add_argument("--enzyme_name","-e",type=str,
                        default="BbsI", required=False,
                        help="""Enzyme name, matching name in the enzyme data file. Default is BbsI""")
    parser.add_argument("--min_oligo_size", "-m", type=int,
                        default=25, required=True,
                        help="""Minimum oligo size required for each DNA sequence. Default 20""")
    parser.add_argument("--max_oligo_size", "-M", type=int,
                        default=100, required=True,
                        help="""Maximum oligo size required for each DNA sequence. Default 100""")
    parser.add_argument("--output", "-o", type=str,
                        default="", required=False,
                        help="""Path to output csv. If not provided, then no csv will be dumped.""")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args=parseargs()

    # Main function
    generate_assembly_library(args.gene_file, args.mutations, args.enzyme_name,
                               args.min_oligo_size, args.max_oligo_size, args.output)
