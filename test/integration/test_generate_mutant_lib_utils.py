import sys
import os
import pytest

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from generate_mutant_lib_utils import (
    make_mut_dict
)

class TestGenerate_Mutant_Lib_Utils:

    @staticmethod
    def test_pass():
        # editable_codon_seq =
        # all_combinations = 
        # name = "name"
        # mutations_df = 

        # make_mut_dict()

        pass

    @staticmethod
    def test_fail():
        pass