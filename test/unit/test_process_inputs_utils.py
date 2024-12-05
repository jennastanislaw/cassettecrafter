import sys
import os
import pytest
from Bio.Seq import Seq
import pandas as pd

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from process_inputs_utils import (
    get_gene_name,
    get_seq,
    biopython_seq_from_str,
    convert_aa_to_dna,
    gen_per_pos_muts,
    split_csl,
    get_allowed_codon_list,
    add_mixed_bases_and_combine,
    get_mixed_base_codon,
    check_leu_arg_ser
)

class TestGet_Gene_Name:
    @staticmethod
    def test_pass():
        type_csv="csv"
        lines_csv=["name1_csv,seq1\n","name2,seq2\n"]
        type_fa="fa"
        lines_fa=[">name1_fa\n","seq line 1\n","seq line 2\n"]
        correct_fa_out="name1_fa"
        correct_csv_out="name1_csv"
        
        fa_out = get_gene_name(lines_fa, type_fa)
        csv_out = get_gene_name(lines_csv, type_csv)

        assert fa_out == correct_fa_out
        assert csv_out == correct_csv_out

    @staticmethod
    def test_fail():
        wrong_type="txt"
        lines_wrong=[[]]

        pytest.raises(ValueError, get_gene_name, lines_wrong, wrong_type)

class TestGet_Seq:
    @staticmethod
    def test_pass():
        pass

    @staticmethod
    def test_fail():
        pass

class TestGen_Per_Pos_Muts:
    @staticmethod
    def test_pass():
        fa_lines=[">genename,Info",
                  "seq1",
                  "seq2"]
        fa_filetype="fa"
        correct_fa_out="seq1seq2"

        csv_lines=["name,fullseq"]
        csv_filetype="csv"
        correct_csv_out="fullseq"
        
        out_fa = get_seq(fa_lines,fa_filetype)
        out_csv = get_seq(csv_lines,csv_filetype)

        assert out_fa == correct_fa_out
        assert out_csv == correct_csv_out

    @staticmethod
    def test_fail():
        lines=["a","b"]
        wrong_filetype = "txt"

        # Should fail because function returns seq and 
        # seq should not be defined if filetype is not allowed
        pytest.raises(UnboundLocalError, get_seq, lines, wrong_filetype)

class TestBiopython_Seq_From_Str:
    @staticmethod

    # Test for DNA. While the function can accept protein (amino acid) input, 
    # this is converted to DNA using convert_aa_to_dna(), and thus that 
    # functionality will not be assessed in this unit test. See convert_aa_to_dna
    # unit tests and integration tests
    def test_pass_dna():
        input_str="ATTGTTCTAGTA"
        input_str_lower=input_str.lower()
        correct_out=Seq("ATTGTTCTAGTA")
        
        output = biopython_seq_from_str(input_str)
        output_lower = biopython_seq_from_str(input_str_lower)

        assert output == correct_out
        assert output_lower == correct_out

    @staticmethod
    def test_fail():
        int_input = 1
        invalid_str="QWERTYUIOP"
        
        pytest.raises(AssertionError, biopython_seq_from_str, int_input)
        pytest.raises(ValueError, biopython_seq_from_str, invalid_str)

class TestConvert_Aa_To_Dna:
    @staticmethod
    def test_pass():
        input_str="APHCVEL"
        input_str_lower="APHCVEL"
        correct_out=Seq("GCACCACACTGCGTAGAACTA")
        
        output = convert_aa_to_dna(input_str)
        output_lower = convert_aa_to_dna(input_str_lower)

        assert output == correct_out
        assert output_lower == correct_out

    @staticmethod
    def test_fail():
        int_input = 1
        invalid_str="QWERTYUIOP"
        
        pytest.raises(AssertionError, convert_aa_to_dna, int_input)
        pytest.raises(ValueError, convert_aa_to_dna, invalid_str)

class TestGen_Per_Pos_Muts:
    @staticmethod
    def test_pass():
        input=pd.DataFrame({"indices":[1,2],
                            "allowed":[["A,C,D"],["F,G"]]})
        correct_out={1:"A,C,D",
                     2:"F,G"}
        
        output= gen_per_pos_muts(input)

        assert output==correct_out
        
    
    @staticmethod
    def test_fail():
        pass

class TestSplit_Csl:
    @staticmethod
    def test_pass():
        input_str="A,B,C"
        correct_output=['A','B','C']

        output=split_csl(input_str)

        assert output==correct_output

'''    @staticmethod
    def test_fail():
        pass

class TestGet_Allowed_Codon_List:
    @staticmethod
    def test_pass():
        pass

    @staticmethod
    def test_fail():
        pass
'''