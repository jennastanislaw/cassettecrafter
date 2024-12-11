import sys
import os
import pandas as pd
import pytest

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from generate_mutant_lib_utils import (
    make_mut_dict
)

class TestGenerate_Mutant_Lib_Utils:

    @staticmethod
    def test_pass_mixed_codon():
        # Mixed codon example
        editable_codon_seq = ['AAG','GAT','GCA','CAA','TCG']
        all_combinations = [("GAT","TAC"),("GAT","CAM"),
                            ("GCA","TAC"),("GCA","CAM")]
        name = "name"

        in_data = {"Position": [2,4],
                "AminoAcids":["A","Y, H"],
                "mut_list":[["A"],["Y","H"]],
                "original":["D", "Q"],
                "codons_original": ["GAT","CAA"],
                "allowed": [["A","D"],["Y","H","Q"]],
                "codons_allowed": [["GAT","GCA"],["TAC","CAM"]]}
        mutations_df = pd.DataFrame.from_dict(in_data, orient="columns").set_index('Position')

        out_data={"name_Q4Y": 'AAGGATGCATACTCG',
                  "name_Q4Q/H": 'AAGGATGCACAMTCG',
                  "name_D2A_Q4Y": 'AAGGCAGCATACTCG',
                  "name_D2A_Q4Q/H": 'AAGGCAGCACAMTCG'}
        
        out = make_mut_dict(editable_codon_seq, all_combinations, name, mutations_df)
        
        assert out == out_data

    @staticmethod
    def test_pass_no_mixed_codon():
        editable_codon_seq = ['AAG','GAT','GCA','CAA','TCG']
        all_combinations = [("GAT","TAC"),("GAT","CCA"), ("GAT","CAA"),
                            ("GCA","TAC"),("GCA","CCA"), ("GCA","CAA")]
        name = "name"

        in_data = {"Position": [2,4],
                "AminoAcids":["A","Y, P"],
                "mut_list":[["A"],["Y","P"]],
                "original":["D", "Q"],
                "codons_original": ["GAT","CAA"],
                "allowed": [["A","D"],["Y","P","Q"]],
                "codons_allowed": [["GAT","GCA"],["TAC","CCA","CAA"]]}
        mutations_df = pd.DataFrame.from_dict(in_data, orient="columns").set_index('Position')

        out_data={"name": 'AAGGATGCACAATCG',
                  "name_Q4Y": 'AAGGATGCATACTCG',
                  "name_Q4P": 'AAGGATGCACCATCG',
                  "name_D2A" :'AAGGCAGCACAATCG',
                  "name_D2A_Q4Y": 'AAGGCAGCATACTCG',
                  "name_D2A_Q4P": 'AAGGCAGCACCATCG'}
        
        out = make_mut_dict(editable_codon_seq, all_combinations, name, mutations_df)
        
        assert out == out_data

    @staticmethod
    def test_fail():
        editable_codon_seq = ['AAG','GAT','GCA','CAA','TCG']
        short_seq = ['AAG','GAT']
        seq_as_str = 'AAGGATGCACAATCG'
        combinations_invalid_codon = [("XXX","TAC")]
        all_combinations = [("GAT","TAC"),("GAT","CCA"), ("GAT","CAA"),
                            ("GCA","TAC"),("GCA","CCA"), ("GCA","CAA")]
        name = "name"
        wrong_name_type = 123

        in_data = {"Position": [2,4],
                "AminoAcids":["A","Y, P"],
                "mut_list":[["A"],["Y","P"]],
                "original":["D", "Q"],
                "codons_original": ["GAT","CAA"],
                "allowed": [["A","D"],["Y","P","Q"]],
                "codons_allowed": [["GAT","GCA"],["TAC","CCA","CAA"]]}
        mutations_df = pd.DataFrame.from_dict(in_data, orient="columns").set_index('Position')

        
        pytest.raises(TypeError, make_mut_dict, editable_codon_seq, 
                      all_combinations, wrong_name_type, mutations_df)
        pytest.raises(TypeError, make_mut_dict, seq_as_str, 
                      all_combinations, name, mutations_df)
        pytest.raises(TypeError, make_mut_dict, editable_codon_seq, 
                      all_combinations, wrong_name_type, mutations_df)
        pytest.raises(KeyError, make_mut_dict, editable_codon_seq, 
                      combinations_invalid_codon, name, mutations_df)
        pytest.raises(IndexError, make_mut_dict, short_seq, 
                      all_combinations, name, mutations_df)
