# loop through data frame of sequences
import pandas as pd
from enzyme_site_replacement import replace_enzyme_site
import enzyme_site_replacement_utils
from process_sequence_list_utils import add_end_restriction_sites, ensure_minimum_length

def replace_enzyme_sites_in_dataframe(df, enzyme_fp, enzyme_name):
    """
    Replace enzyme sites in the DNA sequences located in the 'DNA' column of the DataFrame.

    Parameters:
    - df: DataFrame containing a 'DNA' column.
    - enzyme_fp: File path to the enzyme CSV file.
    - enzyme_name: Name of the enzyme to use for replacement.

    Returns:
    - DataFrame with an additional column 'Modified_DNA' containing the updated sequences.
    """
    # Ensure the 'DNA' column exists
    if 'DNA' not in df.columns:
        raise ValueError("The DataFrame must contain a 'DNA' column.")

    # Define a function to apply to each row
    def replace_site(row):
        DNA = row['DNA']
        try:
            return replace_enzyme_site(enzyme_fp, enzyme_name, DNA)
        except (KeyError, ValueError) as e:
            # Handle exceptions by returning None
            print(f"Error processing DNA: {DNA}. Error: {e}")
            return None  # Does not return sequences that still contain restriction sites

    # Apply the function to each row and create a new column for modified DNA
    df['Modified_DNA'] = df.apply(replace_site, axis=1)

    return df


def process_dna_sequences(df, enzyme_fp, enzyme_name, min_oligo_size):
    """
    Runs the sequence of functions to add restriction sites to each DNA sequence,
    ensures each modified sequence meets the minimum oligo size,
    and generates random oligos if necessary.
    """
    # Step 1: Add enzyme recognition sites
    enzyme_class_objs = enzyme_site_replacement_utils.load_enzymes_from_csv(enzyme_fp)
    enzyme_dict = enzyme_site_replacement_utils.create_enzyme_dict(enzyme_class_objs)

    enzyme = enzyme_dict.get(enzyme_name)

    df = add_end_restriction_sites(df, enzyme_fp, enzyme, min_oligo_size)

    # Step 2: Ensure minimum length for each sequence
    df = ensure_minimum_length(df, min_oligo_size, enzyme)

    return df

if __name__ == "__main__":
    # EXAMPLE USAGE

    # Sample DataFrame with DNA sequences
    data = {
        'DNA': ['ATGCGTACGTAG', 'CGTAGCTAGCAT', 'TGCATGCAATGC']
    }
    df = pd.DataFrame(data)

    # Define the file path to the enzyme CSV file, enzyme name, and minimum oligo size
    enzyme_fp = '../test_data/enzyme_sites.csv'  # Update this path to the location of your enzyme data file
    enzyme_name = 'BbsI'  # Enzyme name as defined in the CSV
    min_oligo_size = 20    # Minimum oligo size required for each DNA sequence

    # Step 1: Replace enzyme sites in each DNA sequence
    df = replace_enzyme_sites_in_dataframe(df, enzyme_fp, enzyme_name)

    # Step 2: Add restriction sites and ensure minimum length
    df = process_dna_sequences(df, enzyme_fp, enzyme_name, min_oligo_size)

    # Display the processed DataFrame
    print(df)