"""
Function to clean up a dataframe input to prepare it to be returned as a 
simplified dataframe. This will also save the cleaned output dataframe to
a csv if a path is provided.

"""
import pandas as pd

def process_output(final_df, output=""):
    """Extract important information from input dataframe, extract columns
        with cassette information, and return the new dataframe. Also writes the
        cleaned dataframe to a csv if an output path is provided.

    Args:
        final_df (Pandas DataFrame): Modified DataFrame with updated DNA sequences
            that have been processed for ordering
        output (str): Path where output csv should be dumped. If empty, no csv produced

    Raises:
        ValueError: If input DataFrame is empty

    Returns:
        Pandas DataFrame: DataFrame with Cassette only columns 
    """
    # Check that input is a datafram 
    if type(final_df) != pd.DataFrame:
        raise ValueError("Input dataframe must be a Pandas DataFrame")
    
    # Assuming final_df is your DataFrame
    final_df['Index'] = range(len(final_df))  # Create an index column with sequential numbers
    final_df.set_index('Index', inplace=True)  # Set the new column as the index

    # Extract only columns starting with "Cassette"
    cassette_columns = [col for col in final_df.columns if col.startswith("Cassette")]

    filtered_df = final_df[cassette_columns]

    # Return error if empty
    if filtered_df.empty:
        raise ValueError(
            "The resulting DataFrame is empty. No columns matching 'Cassette' were found or the data is missing.")
    elif type(output) != str:
        raise TypeError("The path to the output csv must be a string")
    # If path to output file is provided, then dump output csv to that file
    elif output != "":
        filtered_df.to_csv(output,index=False)

    return filtered_df