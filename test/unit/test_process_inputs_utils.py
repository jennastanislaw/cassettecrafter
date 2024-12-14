'''
Unit tests for functions in process_inputs_utils.py, which relate to processesing
and appropriately formatting inputs to the program.

Most functions in process_inputs_utils.py are tested here, with a few notable 
exceptions removed because they make more sense as integration test.
'''

import sys
import os
import pytest
from Bio.Seq import Seq
import pandas as pd

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from process_inputs_utils import (
    # read_gene_input_from_file, # in integrated tests
    get_gene_name,
    get_seq,
    biopython_seq_from_str,
    convert_aa_to_dna,
    # mutation_file_to_df, # in integrated tests
    split_csl,
    #get_allowed_codon_list, #in integrated tests
    add_mixed_bases_and_combine,
    get_mixed_base_codon,
    check_leu_arg_ser
)

class TestGet_Gene_Name:
    """Test get_gene_name(), which should get gene name from the first column of
      csv, or the first line of a fasta
    """
    @staticmethod
    def test_pass():
        """Test csv and fasta inputs. Because this assumes that the lines from the
            file have already been read, the input is a list formatted as it 
            would look when reading from thatf file.
        """
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
        """Ensure correct error is thown when the file type is wrong 
        """
        wrong_type="txt"
        lines=[]

        pytest.raises(ValueError, get_gene_name, lines, wrong_type)

class TestGet_Seq:
    """Tests get_seq() which should get DNA sequence from the second column of 
        csv, or the remaining lines (lines after name) of a fasta
    """
    @staticmethod
    def test_pass():
        """Ensure that function correctly extracts the sequence from a list,
            which represents lines read in from a file. If input is csv, then
            output sequence comes from second coulmn. If it is a fasta then
            output comes from joining all of the lines from line 2 to the end
            of the file
        """
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
        """Ensure correct error is thown when the file type is wrong 
        """
        lines=["a","b"]
        wrong_filetype = "txt"

        # Should fail because function returns seq and 
        # seq should not be defined if filetype is not allowed
        pytest.raises(UnboundLocalError, get_seq, lines, wrong_filetype)

class TestBiopython_Seq_From_Str:
    """Tests biopython_seq_from_str() which converts sequence to a BioPython
        Sequence object. If this is an amnio acid sequence, it will first get
        converted to DNA
    """
    
    # Test for DNA. While the function can accept protein (amino acid) input, 
    # this is converted to DNA using convert_aa_to_dna(), and thus that 
    # functionality will not be assessed in this unit test. See convert_aa_to_dna
    # unit tests 
    @staticmethod
    def test_pass_dna():
        """Ensure that the function correctly accepts a string of DNA and
            appropriately converts it to a BioPython Sequence object with
            the same sequence
        """
        input_str="ATTGTTCTAGTA"
        input_str_lower=input_str.lower()
        correct_out=Seq("ATTGTTCTAGTA")
        
        output = biopython_seq_from_str(input_str)
        output_lower = biopython_seq_from_str(input_str_lower)

        assert output == correct_out
        assert output_lower == correct_out

    @staticmethod
    def test_fail():
        """Ensures that the function throws the correct error when the input
            is not a string, or if the input contains letters that are neither
            amino acids nor nucleotide bases
        """
        int_input = 1
        invalid_str="QWERTYUIOP"
        
        pytest.raises(AssertionError, biopython_seq_from_str, int_input)
        pytest.raises(ValueError, biopython_seq_from_str, invalid_str)

class TestConvert_Aa_To_Dna:
    """Tests function convert_aa_to_dna(), which converts and amino acid string 
        to a DNA sequence composed of codons that correspond to the correct amino acid
        sequence
    """
    @staticmethod
    def test_pass():
        """Ensures that the function correctly converts an uppercase or 
            lowercase amino acid sequence to a DNA sequence
        """
        input_str="APHCVEL"
        input_str_lower="aphcvel"
        correct_out=Seq("GCACCACACTGCGTAGAACTA")
        
        output = convert_aa_to_dna(input_str)
        output_lower = convert_aa_to_dna(input_str_lower)

        assert output == correct_out
        assert output_lower == correct_out

    @staticmethod
    def test_fail():
        """Ensures that the function throws the correct error when the input
            is not a string, or if the input contains letters that are not
            amino acids
        """
        int_input = 1
        invalid_str="QWERTYUIOP"
        
        pytest.raises(AssertionError, convert_aa_to_dna, int_input)
        pytest.raises(ValueError, convert_aa_to_dna, invalid_str)

class TestSplit_Csl:
    """Test function split_csl(), which splits an input string with commas into
        a list of items.
    """
    @staticmethod
    def test_pass():
        """Ensures that the function correclty coverts comma-separated string
            to a list
        """
        input_str="A,B,C"
        correct_output=['A','B','C']

        output=split_csl(input_str)

        assert output==correct_output
    
    @staticmethod
    def test_fail():
        """Ensures that the function throws the correct error when the input
            is not a string
        """
        not_str = 0

        pytest.raises(TypeError, split_csl, not_str)

class TestAdd_Mixed_Bases_And_Combine:
    """Tests function add_mixed_bases_and_combine(), which checks if there are
        any opportunities to introduce mixed base codons based codons in 
        the allowed codon list which have the same first two bases
    """
    @staticmethod
    def test_pass():
        """Ensures that function correctly identifies opportunities for mixed
            base codon to be used and updated the allowed codon list accordingly
        """
        input_codons_one = ["ATG", "ATC"]
        expected_output_one = ["ATS"]

        input_codons_two = ["ATG", "ATC", "CCC","GAT","GAA"]
        expected_output_two = ["ATS", "CCC", "GAW"]

        no_change_input=["ATC", "CCC"]
        expected_output_no_change = ["ATC", "CCC"]

        empty=[]
        expected_empty_out=[]

        output_one = add_mixed_bases_and_combine(input_codons_one)
        output_two = add_mixed_bases_and_combine(input_codons_two)
        output_no_change = add_mixed_bases_and_combine(no_change_input)
        output_empty = add_mixed_bases_and_combine(empty)

        assert expected_output_one.sort() == output_one.sort()
        assert expected_output_two.sort() == output_two.sort()
        assert output_empty == expected_empty_out
        assert expected_output_no_change == no_change_input

    @staticmethod
    def test_fail():
        """Ensures that the function throws the correct error when the input
            is not a list
        """
        wrong_type=0

        pytest.raises(TypeError, add_mixed_bases_and_combine, wrong_type)

class TestGet_Mixed_Base_Codon:
    """Tests get_mixed_base_codon(), which find the appropriate codon containing
        a mixed base in the third position which satisfies both of the input codons. 
    """
    @staticmethod
    def test_pass():
        """Ensures that the correct mixed base codon is returned
        """
        input1 = "ATT"
        input2 = "ATG"
        correct_output="ATK"

        output = get_mixed_base_codon(input1, input2)

        assert output == correct_output

    @staticmethod
    def test_fail():
        """Ensures that the correct error is thrown if input codon is the
            wrong length, type, or the first two bases of the input codons do not
            match 
        """
        wrong_type=0
        wrong_len="ATAGCTGA"
        normal_codon="ATG"
        no_match="CAA"

        pytest.raises(TypeError, get_mixed_base_codon, wrong_type, normal_codon)
        pytest.raises(ValueError, get_mixed_base_codon, normal_codon, wrong_len)
        pytest.raises(ValueError, get_mixed_base_codon, normal_codon, no_match)


class TestCheck_Leu_Arg_Ser:
    """Tests check_leu_arg_ser(), which catches and handles edge cases when looking
        for opportunities to introduce mixed base codons. See the function docstring
        for more specific details about what this means
    """
    @staticmethod
    def test_pass():
        """Ensures that each of these edge case AAs are appropriately handled,
            and that there is no change made if the input base is not one of the
            exceptions.
        """
        leu_input = "CTT"
        arg_input = "CGT"
        ser_input = "TCT"
        none_input = "ATT"

        correct_leu_output = "TTA"
        correct_arg_output = "AGA"
        correct_ser_output = "AGC"
        correct_none_output = "ATT"

        out_leu = check_leu_arg_ser(leu_input)
        out_arg = check_leu_arg_ser(arg_input)
        out_ser = check_leu_arg_ser(ser_input)
        out_none = check_leu_arg_ser(none_input)
        
        assert correct_leu_output == out_leu
        assert correct_arg_output == out_arg
        assert correct_ser_output == out_ser
        assert correct_none_output == out_none

    @staticmethod
    def test_fail():
        """Ensures that the correct error is thrown if input codon is the
            wrong length or type
        """
        wrong_type=0
        wrong_len="ATAGCTGA"

        pytest.raises(ValueError, check_leu_arg_ser, wrong_len)
        pytest.raises(TypeError, check_leu_arg_ser, wrong_type)
