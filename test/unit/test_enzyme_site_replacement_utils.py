import sys
import os

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

# Now you can import your modules
from enzyme_site_replacement_utils import (
    Enzyme,
    create_enzyme_dict,
    load_enzymes_from_csv,
    find_matching_sites,
    generate_synonymous_codons_dna,
    get_codons_with_full_indices,
    get_affected_codons_by_recognition_sites,
)

import pytest

# Sample data for testing
@pytest.fixture
def enzyme_data():
    return [
        Enzyme(name='EcoRI', fwd_recognition_site='GAATTC', rev_recognition_site='CTTAAG', spacer_length=2, OH_length=4),
        Enzyme(name='BamHI', fwd_recognition_site='GGATCC', rev_recognition_site='CCTAGG', spacer_length=2, OH_length=4)
    ]

def test_create_enzyme_dict(enzyme_data):
    enzyme_dict = create_enzyme_dict(enzyme_data)
    assert len(enzyme_dict) == 2
    assert enzyme_dict['EcoRI'].fwd_recognition_site == 'GAATTC'
    assert enzyme_dict['BamHI'].rev_recognition_site == 'CCTAGG'

def test_load_enzymes_from_csv(tmpdir):
    # Create a temporary CSV file for testing
    csv_content = """Enzyme,Fwd_recognition_site,Rev_recognition_site,spacer_length,OH_length
EcoRI,GAATTC,CTTAAG,0,4
BamHI,GGATCC,CCTAGG,0,4
"""
    csv_file = tmpdir.join("enzymes.csv")
    with open(csv_file, 'w') as f:
        f.write(csv_content)

    enzymes = load_enzymes_from_csv(str(csv_file))
    assert len(enzymes) == 2
    assert enzymes[0].name == 'EcoRI'
    assert enzymes[1].fwd_recognition_site == 'GGATCC'

def test_find_matching_sites():
    enzyme = Enzyme(name='EcoRI', fwd_recognition_site='GAATTC', rev_recognition_site='CTTAAG', spacer_length=0, OH_length=4)
    dna_sequence = "AGCTGAATTCAGCTCTTAAGCTAGATG"

    fwd_matches, rev_matches = find_matching_sites(enzyme, dna_sequence)
    assert fwd_matches == [4]
    assert rev_matches == [14]

class TestCodonFunctions:
    def test_generate_synonymous_codons_dna_valid(self):
        """Test for valid synonymous codons."""
        codon = 'TTT'  # For Phenylalanine (F)
        synonymous_codons = generate_synonymous_codons_dna(codon)

        assert 'TTC' in synonymous_codons  # Check if 'TTC' is a synonymous codon
        assert len(synonymous_codons) == 1  # Verify the number of synonymous codons

    def test_generate_synonymous_codons_dna_invalid(self):
        """Test for handling of invalid codons."""
        with pytest.raises(ValueError):
            generate_synonymous_codons_dna('INVALID')  # Expect a ValueError for invalid codon

def test_get_codons_with_full_indices():
    dna_sequence = "ATGTTTGCAATG"
    codons_with_indices = get_codons_with_full_indices(dna_sequence)

    assert len(codons_with_indices) == 4  # 4 complete codons
    assert codons_with_indices[0] == ('ATG', [0, 1, 2])
    assert codons_with_indices[1] == ('TTT', [3, 4, 5])
    assert codons_with_indices[2] == ('GCA', [6, 7, 8])
    assert codons_with_indices[3] == ('ATG', [9, 10, 11])

def test_get_affected_codons_by_recognition_sites():
    dna_sequence = "ATGGAATTCGCTCTTAAGGAT"
    recognition_idxs = [3, 14]
    recognition_site = 'GAATTC'

    affected_codons = get_affected_codons_by_recognition_sites(dna_sequence, recognition_idxs, recognition_site)

    assert len(affected_codons) == 2  # There should be affected codons for both sites
    assert affected_codons[3] == [('GAA', [3, 4, 5]), ('TTC', [6, 7, 8])]  # Site at index 3
    assert affected_codons[14] == [('CTT', [12, 13, 14]), ('AAG', [15, 16, 17]), ('GAT', [18, 19, 20])]  # Site at index 14