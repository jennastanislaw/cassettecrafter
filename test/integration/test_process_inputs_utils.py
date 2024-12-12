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
    @staticmethod
    def test_pass(request):
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
        wrong_type = 1
        wrong_path="doesnt_exist.txt"

        pytest.raises(TypeError, read_gene_input_from_file, wrong_type)
        pytest.raises(ValueError, read_gene_input_from_file, wrong_path)

class TestMutation_File_To_Df:
    @staticmethod
    def test_pass(request):
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
        wrong_type = 1
        wrong_extension = "file.pdb"
        seq=Seq("ATC")
        wrong_seq="ATC"

        pytest.raises(TypeError, mutation_file_to_df, wrong_type, seq)
        pytest.raises(TypeError, mutation_file_to_df, wrong_extension, seq)
        pytest.raises(TypeError, mutation_file_to_df, "file.csv", wrong_seq)
    
class TestGet_Allowed_Codon_List:
    @staticmethod
    def test_pass():
        mutations_input = [["Y","H"]]
        codons_original = ["CAA"]

        expected_output = [["TAC","CAM"]]

        output = get_allowed_codon_list(mutations_input,codons_original)

        assert output == expected_output

    @staticmethod
    def test_fail():
        non_codon_input = [["XFT"]]
        codons_original = ["CAA"]

        mutation_input_long=[[None],[None]]
        notlist=1

        pytest.raises(KeyError, get_allowed_codon_list, non_codon_input,codons_original)
        pytest.raises(TypeError,get_allowed_codon_list,non_codon_input,notlist)
        pytest.raises(ValueError, get_allowed_codon_list, mutation_input_long,codons_original)