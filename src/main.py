import argparse
# import utils
from helper_functions import *

"""
Currently running with:
python3 main.py -m ../test_data/demo_mutation_list.csv -f ../test_data/LY011_test_seq_single.csv
python3 main.py -m ../test_data/demo_mutation_list.csv -f ../test_data/LY011_test_seq_single.fa

"""

def generate_assembly(gene_file, mutations, backbone):
    name, starting_dna = read_input(gene_file)
    # Do some checks on the input in above fn
   
    mutations_df = mutation_file_to_df(mutations, starting_dna)
    print(mutations_df)
    print(starting_dna)

    library_dict = generate_mutant_lib(starting_dna,mutations_df, name)
    print(library_dict)
    #generate_dna_lib

def parseargs():
    #TODO: add docstring and decriptions 
    parser=argparse.ArgumentParser(
        description="""""", 
        prog="parse_inputs.py"
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

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args=parseargs()
    # Main function
    generate_assembly(args.gene_file, args.mutations, args.backbone)


