import sys
import os
import pytest

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from generate_mutant_lib_utils import (
    # make_mut_dict, # in intergation tests
    split_to_codons
)

class TestSplit_To_Codons:
    @staticmethod
    def test_pass():
        input="ABCDEFGHI"
        expected_output=["ABC","DEF","GHI"]

        output = split_to_codons(input)

        assert output == expected_output

    @staticmethod
    def test_fail():
        not_correct_len="ACGT"

        pytest.raises(ValueError, split_to_codons, not_correct_len)
