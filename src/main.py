"""
From provided sequence and allowed mutations, produces a library of sequences
that contain all combinations of the desired mutations that are compatible with 
Golden Gate Assembly

Usage: python3 main.py -m [allowed mutations file] -b [plasmid backbone] 
        -e [restriction enzyme name] -d [file with enzyme data] 
        -m [minimum oligo size] -f [file with genes to insert]

Example: python3 main.py -m ./test_data/demo_mutation_list.csv 
                -f ./test_data/LY011_test_seq_single.csv
    When run from the directory above src. Adjust paths as needed

"""

import argparse
import os

from process_inputs import process_inputs
from generate_mutant_lib import generate_mutant_lib
from process_sequence_list import (replace_enzyme_sites_in_dataframe, 
                                    process_dna_sequences)

def generate_assembly_library(gene_file, mutations, backbone, enzyme_data, 
                      enzyme_name, min_oligo_size):

    name, starting_dna, mutations_df = process_inputs(gene_file,mutations)

    library_df = generate_mutant_lib(starting_dna,mutations_df, name)

    replace_enzyme_sites_in_dataframe(library_df, enzyme_data, enzyme_name)

    final_df = process_dna_sequences(library_df, enzyme_data, enzyme_name,
                                     min_oligo_size)

    # print(final_df)

    # assert final_df.equals(library_df)

    return final_df

def parseargs():
    #TODO: add docstring and decriptions 
    parser=argparse.ArgumentParser(
        description="""""", 
        prog="generate_assembly_library.py"
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