import sys
import os
import pytest

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from process_inputs_utils import (
    get_gene_name,
    get_dna_seq,
    gen_per_pos_muts,
    split_csl,
    get_allowed_codon_list,
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

'''class TestGet_Dna_Seq:
    @staticmethod
    def test_pass():
        pass

    @staticmethod
    def test_fail():
        pass

class TestGen_Per_Pos_Muts:
    @staticmethod
    def test_pass():
        pass

    @staticmethod
    def test_fail():
        pass
'''
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