import sys
import os
import pytest
import pandas as pd

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from process_output import (
    process_output
)

class TestProcess_Output:

    @staticmethod
    def test_pass(tmp_path):
        df=pd.DataFrame([[1,"ABC","abc",0],[2,"DEF","def",0]],
                        columns=["name","Cassette1","Cassette2","extra"])
        corr_out=pd.DataFrame([["ABC","abc"],["DEF","def"]],
                    columns=["Cassette1","Cassette2"])
        out_csv = "temp.csv"
        temp_file=f"{tmp_path}/{out_csv}"

        output = process_output(df)
        process_output(df,temp_file)

        assert corr_out.equals(output)

        assert os.path.exists(temp_file)
        assert corr_out.equals(pd.read_csv(temp_file,index_col=False))

    @staticmethod
    def test_fail():
        not_df="not_df"
        empty_df=pd.DataFrame([])
        real_df=pd.DataFrame([[1,2,3],["a","b","c"]],
                             columns=["name","Cassette1","Cassette2"])
        output_int=1

        pytest.raises(ValueError, process_output, not_df)
        pytest.raises(ValueError, process_output, empty_df)
        pytest.raises(TypeError, process_output, real_df, output_int)