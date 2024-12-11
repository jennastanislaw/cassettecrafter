"""
This script processes a DataFrame of DNA sequences by adding terminal restriction sites,
ensuring sequences meet a specified minimum oligo size, and modifying sequences as needed.

Dependencies:
    pandas (pd): Required for handling the DataFrame containing DNA sequences.
    process_sequence_list_utils (add_end_restriction_sites, ensure_minimum_length):
        - add_end_restriction_sites: Adds enzyme recognition sites to DNA sequences.
        - ensure_minimum_length: Ensures sequences meet the specified minimum oligo size.

Functions:
    process_dna_sequences(df, enzyme, min_oligo_size):
        Processes DNA sequences by applying terminal restriction sites, checking
        for minimum size compliance, and making necessary adjustments. Returns
        a modified DataFrame.
"""

from utils.process_sequence_list_utils import add_end_restriction_sites, ensure_minimum_length


def process_dna_sequences(df, enzyme, min_oligo_size):
    """
    Processes a DataFrame of DNA sequences by adding terminal restriction sites,
    ensuring minimum oligo size, and modifying sequences as necessary.

    Args:
        df (pd.DataFrame): Input DataFrame containing DNA sequences.
        enzyme (Enzyme): Enzyme object with recognition site properties.
        min_oligo_size (int): Minimum allowed size for oligo sequences.

    Returns:
        pd.DataFrame: Modified DataFrame with updated DNA sequences.
    """
    # Identify columns containing cassette DNA sequences
    cassette_columns = [col for col in df.columns if col.startswith('Cassette')]

    # Iterate through cassette columns to process sequences
    for col in cassette_columns:
        # Add restriction sites to the sequences in the current column
        df = add_end_restriction_sites(df, col, enzyme)

        # Ensure sequences meet the minimum oligo size
        df = ensure_minimum_length(df, col, min_oligo_size, enzyme)

    return df

