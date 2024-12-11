import sys
import os
import pandas as pd
from Bio.Seq import Seq

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from generate_mutant_lib import (
    generate_mutant_lib
)

class TestGenerate_Mutant_Lib:
    @staticmethod
    def test_pass():
        og_codon_seq = Seq("AAGGATGCACAATCG")

        in_data = {"Position": [2,4],
                "AminoAcids":["A","Y, P"],
                "mut_list":[["A"],["Y","P"]],
                "original":["D", "Q"],
                "codons_original": ["GAT","CAA"],
                "allowed": [["A","D"],["Y","P","Q"]],
                "codons_allowed": [["GAT","GCA"],["TAC","CCA","CAA"]]}
        mutations_df = pd.DataFrame.from_dict(in_data, orient="columns").set_index('Position')

        name="my_gene"

        expected_out=pd.DataFrame({"name":["my_gene", "my_gene_Q4Y",
                                           "my_gene_Q4P", "my_gene_D2A",
                                           "my_gene_D2A_Q4Y", "my_gene_D2A_Q4P"],
                                   "DNA": ['AAGGATGCACAATCG', 'AAGGATGCATACTCG',
                                           'AAGGATGCACCATCG','AAGGCAGCACAATCG',
                                           'AAGGCAGCATACTCG','AAGGCAGCACCATCG']}) 

        out = generate_mutant_lib(og_codon_seq,mutations_df, name)

        out_sorted = out.sort_values(by=['name']).reset_index(drop=True)
        expected_out_sorted = expected_out.sort_values(by=['name']).reset_index(drop=True)

        assert (out_sorted.equals(expected_out_sorted))
