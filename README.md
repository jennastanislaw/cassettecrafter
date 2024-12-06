# cassettecrafter
This software generates a library of DNA sequences casettes that contains
all possible combinations of allowed mutations at specified positions, with the necessary recognition sequences for Golden Gate Assembly. Note that the GUI portion of this code and the DNA processing portion of the code are not yet integrated. Please see the GUI directory for its own README and instructions for use.

## Usage
From the src directory, sequences can be generated by running main.py

Usage: `python3 main.py -m [allowed mutations file] -b [plasmid backbone] -e [restriction enzyme name] -d [file with enzyme data] -m [minimum oligo size] -f [file with genes to insert]`

To run an example from the `cassettecrafter` directory (parent directory of `src`):`python3 main.py -m ./test_data/demo_mutation_list.csv -f ./test_data/LY011_test_seq_single.csv`

### Proper formatting of input:
* gene_file : Gene sequence (currently a codon sequence) can be provided as a fasta file in the traditional format or a csv without a header, with each row formatted as [gene_name],[DNA_sequence]
    * cassettecrafter/test_data/LY011_test_seq_single.fa
    * cassettecrafter/test_data/LY011_test_seq_single.csv
* mutations : Allowed mutations can currently only be provided as a csv, with a header row. The first column should correspond to the position in the amino acid sequence (1-indexed), and the second column should contain a string of allowed amino acids mutations at that positions, separated by commas.
    * cassettecrafter/test_data/demo_mutation_list.csv
* enzyme_data : Information about restriction enzymes used for Golden Gate Assembly. Can contain data about multiple enzymes, but only one will be used for inserting recognition and cut sites into the provided gene, based on the user-provided enzyme name. Should be provided as a csv with a header row, and the following columns: Enzyme,Fwd_recognition_site,Rev_recognition_site,spacer_length,OH_length
    * cassettecrafter/test_data/enzyme_sites.csv

## Installation and Setup
After cloning the GitHub repository, the user can create the conda environment containing packages needed to run all of the sofware with: 
 `conda env create -f environment.yml`

The following Python packages are required for this software and are included:
* Pandas
* BioPython
* Optional: pytest


## Notes for troubleshooting
Indexing for amino acid mutations starts at 1, not at 0 