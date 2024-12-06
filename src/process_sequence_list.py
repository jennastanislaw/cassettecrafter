# loop through data frame of sequences
from process_sequence_list_utils import add_end_restriction_sites, ensure_minimum_length


def process_dna_sequences(df, enzyme, min_oligo_size):
    """
    Runs the sequence of functions to add restriction sites to each DNA sequence,
    ensures each modified sequence meets the minimum oligo size,
    and generates random oligos if necessary.
    """
    # Step 1: Add enzyme recognition sites
    cassette_columns = [col for col in df.columns if col.startswith('Cassette')]

    for col in cassette_columns:
        df = add_end_restriction_sites(df, col, enzyme)

        # Step 2: Ensure minimum length for each sequence
        df = ensure_minimum_length(df, col, min_oligo_size, enzyme)

    return df

