""" 
Integration tests for process_inputs.py and its main function. This modules serves
to take an input sequence and mutation csv and return the correct seuqnce name,
 sequence, and mutation dataframe.
"""

import sys
import os
from Bio.Seq import Seq
import pandas as pd

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from process_inputs import (
    process_inputs
)

class TestProcess_Inputs:
    """Test process_inputs(), which takes an input sequence and mutation csv and
        returns the correct seuqnce name, sequence, and mutation dataframe.
        Passing functions here also helps verify that helper functions
        read_gene_input_from_file, biopython_seq_from_str, 
        and mutation_file_to_df are working correctly
    """

    # Define mutation data. This is what the dataframe should look like after
    # parsing input mutation csv and input sequence
    mut_data = {"Position": [4],
                "AminoAcids":["Y, H"],
                "mut_list":[["Y","H"]],
                "original":["L"],
                "codons_original": ["TTA"],
                "allowed": [["Y","H","L"]],
                "codons_allowed": [["TAC","CAC","TTA"]]}
        
    mut_out = pd.DataFrame.from_dict(mut_data, orient="columns").set_index('Position')

    @staticmethod
    def test_from_csv(request):
        """Ensures that the function can take a csv sequence input and mutation
            csv and produce the correct output gene name, sequence, and mutation 
            datafame. 

         Args:
            request : pytest fixture used to get the path to the csv used for 
            testing, relative to the root directory from which pytest is run
        """
        rootdir=request.config.rootdir
        csv_input= f"{rootdir}/test_data/test_seq_single.csv"
        mut_csv = f"{rootdir}/test/integration/sample_mut.csv"

        expected_out=("example_1",
                      Seq("AAGGATGCATTATCGAAAGCCGTGAGTGAACGGAGAGGTCAAGTTGGCGGCTGAAAGTCGATTTATTGACGGCAGAAATTGACTCTGTCAGTCCTTTCGTAGTTGTGGCGACTTCAACATGGTACA"),
                      TestProcess_Inputs.mut_out)
        
        out = process_inputs(csv_input,mut_csv)

        assert expected_out[0] == out[0]
        assert expected_out[1] == out[1]
        assert expected_out[2].equals(out[2])

    @staticmethod
    def test_from_fa(request):
        """Ensures that the function can take a fasta sequence input and mutation
            csv and produce the correct output gene name, sequence, and mutation 
            datafame. 

         Args:
            request : pytest fixture used to get the path to the fasta used for 
            testing, relative to the root directory from which pytest is run
        """
        rootdir=request.config.rootdir
        fa_input=f"{rootdir}/test_data/test_seq_single.fa"
        mut_csv = f"{rootdir}/test/integration/sample_mut.csv"

        expected_out=("example_1",
                      Seq("AAGGATGCATTATCGAAAGCCGTGAGTGAACGGAGAGGTCAAGTTGGCGGCTGAAAGTCGATTTATTGACGGCAGAAATTGACTCTGTCAGTCCTTTCGTAGTTGTGGCGACTTCAACATGGTACA"),
                      TestProcess_Inputs.mut_out)
        
        out = process_inputs(fa_input,mut_csv)

        assert expected_out[0] == out[0]
        assert expected_out[1] == out[1]
        assert expected_out[2].equals(out[2])

    @staticmethod
    def test_from_string(request):
        """Ensures that the function can take a string as sequence input and mutation
            csv and produce the correct output gene name, sequence, and mutation 
            datafame. 
        """
        rootdir=request.config.rootdir
        seq_str="AAGGATGCATTATCGAAAGCCGTGAGTGAACGGAGAGGTCAAGTTGGCGGCTGAAAGTCGATTTATTGACGGCAGAAATTGACTCTGTCAGTCCTTTCGTAGTTGTGGCGACTTCAACATGGTACA"
        mut_csv = f"{rootdir}/test/integration/sample_mut.csv"

        expected_out=("gene1", # this is hard-coded into the function
                      Seq("AAGGATGCATTATCGAAAGCCGTGAGTGAACGGAGAGGTCAAGTTGGCGGCTGAAAGTCGATTTATTGACGGCAGAAATTGACTCTGTCAGTCCTTTCGTAGTTGTGGCGACTTCAACATGGTACA"),
                      TestProcess_Inputs.mut_out)
        
        out = process_inputs(seq_str,mut_csv)

        assert expected_out[0] == out[0]
        assert expected_out[1] == out[1]
        assert expected_out[2].equals(out[2])