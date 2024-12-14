
""" 
Integration tests for process_inputs_utils.py, specifically the functions 
read_gene_input_from_file(),  mutation_file_to_df(), get_allowed_codon_list(). 
"""

import sys
import os
import pandas as pd
from Bio.Seq import Seq
import pytest

# Add the src and src/util directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from process_inputs_utils import (
    read_gene_input_from_file, 
    mutation_file_to_df, 
    get_allowed_codon_list, 
)

class TestRead_Gene_Input_From_File:
    """Tests read_gene_input_from_file(), which extract gene name and sequence 
        from an input file
    """
    @staticmethod
    def test_pass(request):
        """Ensures that the csv or fasta input is parsed correctly and produces a
            tuple with the name and sequence from the file. Passing test helps
            verify that helper functions get_gene_name and get_seq are working
            correctly
        Args:
            request : pytest fixture used to get the path to the csv or fasta used for 
            testing, relative to the root directory from which pytest is run
        """
        rootdir=request.config.rootdir
        csv_input= f"{rootdir}/test_data/test_seq_single.csv"
        fa_input=f"{rootdir}/test_data/test_seq_single.fa"

        expected_out= ("example_1",
            """AAGGATGCATTATCGAAAGCCGTGAGTGAACGGAGAGGTCAAGTTGGCGGCTGAAAGTCGATTTATTGACGGCAGAAATTGACTCTGTCAGTCCTTTCGTAGTTGTGGCGACTTCAACATGGTACA""")
        
        out_csv = read_gene_input_from_file(csv_input)
        out_fa = read_gene_input_from_file(fa_input)

        assert out_csv == expected_out
        assert out_fa == expected_out

    @staticmethod
    def test_fail():
        """Ensures that the correct errors are thrown when the input is the wrong
            type or points to a path that doesn't exist
        """
        wrong_type = 1
        wrong_path="doesnt_exist.txt"

        pytest.raises(TypeError, read_gene_input_from_file, wrong_type)
        pytest.raises(ValueError, read_gene_input_from_file, wrong_path)

class TestMutation_File_To_Df:
    """Tests mutation_file_to_df(), which accepts csv input containing columns of
        positions in an amino acid sequence and what they can be converted to,
        and returns a dataframe which a variety of mutation information, based on
        the mutation file and the provided sequence
    """
    @staticmethod
    def test_pass(request):
        """Ensures that the expected dataframe is produced. Passing test helps
            verify that the helper function get_allowed_codon_list is working
            correctly
        Args:
            request : pytest fixture used to get the path to the csv used for 
            testing, relative to the root directory from which pytest is run
        """
        rootdir=request.config.rootdir
        mut_csv = f"{rootdir}/test/integration/sample_mut.csv"
        input_seq=Seq("ATGAATCTACAA")
        data = {"Position": [4],
                "AminoAcids":["Y, H"],
                "mut_list":[["Y","H"]],
                "original":["Q"],
                "codons_original": ["CAA"],
                "allowed": [["Y","H","Q"]],
                "codons_allowed": [["TAC","CAM"]]}
        
        expected_out = pd.DataFrame.from_dict(data, orient="columns").set_index('Position')

        out = mutation_file_to_df(mut_csv,input_seq)

        assert (out.equals(expected_out))

    @staticmethod
    def test_fail():
        """Ensures that the correct error type is thrown when inputs are not the
            correct type
        """
        wrong_type = 1
        wrong_extension = "file.pdb"
        seq=Seq("ATC")
        wrong_seq="ATC"

        pytest.raises(TypeError, mutation_file_to_df, wrong_type, seq)
        pytest.raises(TypeError, mutation_file_to_df, wrong_extension, seq)
        pytest.raises(TypeError, mutation_file_to_df, "file.csv", wrong_seq)
    
class TestGet_Allowed_Codon_List:
    """Tests get_allowed_codon_list(), which converts amino acid information 
        to codon information for each of the allowed mutations
    """
    @staticmethod
    def test_pass():
        """Ensures that the correct list of allowed codons is returned from the 
            allowed amino acid input. This includes taking into account the
            need to add mixed base codons if necessary. Passing test helps
            verify that helper function add_mixed_bases_and_combine is working 
            correctly.
        """
        mutations_input = [["Y","H"]]
        codons_original = ["CAA"]

        expected_output = [["TAC","CAM"]]

        output = get_allowed_codon_list(mutations_input,codons_original)

        assert output == expected_output

    @staticmethod
    def test_fail():
        """Ensures that correct errors are thrown when inputs are not the correct
            type, the list of allowed amino acid mutations has a letter that 
            doesn't correspond to a valid amino acid, or the mutation list is 
            shorter than sequence to be mutated
        """
        non_aa_input = [["XFT"]]
        codons_original = ["CAA"]

        mutation_input_long=[[None],[None]]
        notlist=1

        pytest.raises(KeyError, get_allowed_codon_list, non_aa_input,codons_original)
        pytest.raises(TypeError,get_allowed_codon_list,non_aa_input,notlist)
        pytest.raises(ValueError, get_allowed_codon_list, mutation_input_long,codons_original)