""" 
Unit tests for generate_mutant_lib_utils.py, which help process information from
the mutant dataframe.

"""

import sys
import os
import pytest
import pandas as pd

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from generate_mutant_lib_utils import (
    make_mut_dict,
    split_to_codons
)

class TestMake_Mut_Dict:
    """Tests make_mut_dict(), which take mutation information and sequence info and
        generates a dictionary where the first coulmn is the mutation name and
        the second column is the sequence of that mutant. Passing tests also help
        verify that the helper functio 
    """

    @staticmethod
    def test_pass_mixed_codon():
        """Ensures that all of the correct names and sequences are generated when a
            mixed base codon is included in the list of allowed codons. This is a
            special case because the name of both of the amino acids that the 
            mixed base codon is encoding must be extracted from the mixed base 
            codon alone.
        """
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
        """Ensures that all of the correct names and sequences are generated 
            from the input seuqence and the allowed mutation information
        """
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
        """Ensures that the correct error is thrown when the variables are the 
            wrong type, a codon in the allowed codons is not valid, or a position
            where a mutation is allow is not within the length of the sequence
            that is being mutated
        """
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


class TestSplit_To_Codons:
    """Tests split_to_codons(), which splits a DNA sequence to a list of codons
        by dividing it into fragments of length 3
    """
    @staticmethod
    def test_pass():
        """Ensures that input is properly split into sets of three 
        """
        input="ABCDEFGHI"
        expected_output=["ABC","DEF","GHI"]

        output = split_to_codons(input)

        assert output == expected_output

    @staticmethod
    def test_fail():
        """Ensures that correct error is thrown if input is not a string or
            if its length is not a multiple of 3
        """
        not_correct_len="ACGT"
        wrong_type=0

        pytest.raises(ValueError, split_to_codons, not_correct_len)
        pytest.raises(TypeError, split_to_codons, wrong_type)
