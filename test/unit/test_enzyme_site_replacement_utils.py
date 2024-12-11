"""
This script contains tests for functions that handle DNA sequence manipulation
with respect to enzyme recognition sites and codon generation. The tests cover
enzyme loading, enzyme site matching, codon-related functions, and the identification
of affected codons by enzyme recognition sites.

Functions:
- enzyme_data: Provides sample data for testing enzyme functions.
- test_create_enzyme_dict: Tests the creation of a dictionary of enzyme objects.
- test_load_enzymes_from_csv: Tests loading enzymes from a CSV file.
- test_find_matching_sites: Tests finding matching enzyme recognition sites in a DNA sequence.
- test_generate_synonymous_codons_dna_valid: Tests generating valid synonymous codons for a given codon.
- test_generate_synonymous_codons_dna_invalid: Tests handling of invalid codons.
- test_get_codons_with_full_indices: Tests extracting codons and their indices from a DNA sequence.
- test_get_affected_codons_by_recognition_sites: Tests identifying affected codons by enzyme recognition sites.

Dependencies:
- `sys`: To modify the Python path and include the src directory.
- `os`: For path manipulation.
- `enzyme_site_replacement_utils`: Contains functions for enzyme site recognition and codon operations.
- `pytest`: For running tests and assertions.
"""

import sys
import os
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/utils')))

from enzyme_site_replacement_utils import (
    Enzyme,
    create_enzyme_dict,
    load_enzymes_from_csv,
    find_matching_sites,
    generate_synonymous_codons_dna,
    get_codons_with_full_indices,
    get_affected_codons_by_recognition_sites,
)


@pytest.fixture
def enzyme_data():
    """Fixture that provides sample enzyme data for testing."""
    return [
        Enzyme(name='EcoRI', fwd_recognition_site='GAATTC', rev_recognition_site='CTTAAG', spacer_length=2,
               OH_length=4),
        Enzyme(name='BamHI', fwd_recognition_site='GGATCC', rev_recognition_site='CCTAGG', spacer_length=2, OH_length=4)
    ]


def test_create_enzyme_dict(enzyme_data):
    """Test the creation of an enzyme dictionary."""
    enzyme_dict = create_enzyme_dict(enzyme_data)
    # Check if the enzyme dictionary is correctly populated
    assert len(enzyme_dict) == 2
    assert enzyme_dict['EcoRI'].fwd_recognition_site == 'GAATTC'
    assert enzyme_dict['BamHI'].rev_recognition_site == 'CCTAGG'


def test_load_enzymes_from_csv(tmpdir):
    """Test loading enzymes from a CSV file."""
    # Create a temporary CSV file for testing
    csv_content = """Enzyme,Fwd_recognition_site,Rev_recognition_site,spacer_length,OH_length
EcoRI,GAATTC,CTTAAG,0,4
BamHI,GGATCC,CCTAGG,0,4
"""
    csv_file = tmpdir.join("enzymes.csv")
    with open(csv_file, 'w') as f:
        f.write(csv_content)

    enzymes = load_enzymes_from_csv(str(csv_file))
    # Ensure the enzymes are correctly loaded from CSV
    assert len(enzymes) == 2
    assert enzymes[0].name == 'EcoRI'
    assert enzymes[1].fwd_recognition_site == 'GGATCC'


def test_find_matching_sites():
    """Test finding matching forward and reverse recognition sites in a DNA sequence."""
    enzyme = Enzyme(name='EcoRI', fwd_recognition_site='GAATTC', rev_recognition_site='CTTAAG', spacer_length=0,
                    OH_length=4)
    dna_sequence = "AGCTGAATTCAGCTCTTAAGCTAGATG"

    fwd_matches, rev_matches = find_matching_sites(enzyme, dna_sequence)
    # Assert correct matches for both forward and reverse sites
    assert fwd_matches == [4]
    assert rev_matches == [14]


class TestCodonFunctions:
    """Tests for codon-related functions."""

    def test_generate_synonymous_codons_dna_valid(self):
        """Test for valid synonymous codon generation."""
        codon = 'TTT'  # For Phenylalanine (F)
        synonymous_codons = generate_synonymous_codons_dna(codon)

        assert 'TTC' in synonymous_codons  # Check if 'TTC' is a synonymous codon
        assert len(synonymous_codons) == 1  # Verify the number of synonymous codons

    def test_generate_synonymous_codons_dna_invalid(self):
        """Test handling of invalid codons."""
        with pytest.raises(ValueError):
            generate_synonymous_codons_dna('INVALID')  # Expect a ValueError for invalid codon


def test_get_codons_with_full_indices():
    """Test extracting codons and their indices from a DNA sequence."""
    dna_sequence = "ATGTTTGCAATG"
    codons_with_indices = get_codons_with_full_indices(dna_sequence)

    # Check that the correct number of codons and their indices are returned
    assert len(codons_with_indices) == 4  # 4 complete codons
    assert codons_with_indices[0] == ('ATG', [0, 1, 2])
    assert codons_with_indices[1] == ('TTT', [3, 4, 5])
    assert codons_with_indices[2] == ('GCA', [6, 7, 8])
    assert codons_with_indices[3] == ('ATG', [9, 10, 11])


def test_get_affected_codons_by_recognition_sites():
    """Test extracting codons affected by enzyme recognition sites."""
    dna_sequence = "ATGGAATTCGCTCTTAAGGAT"
    recognition_idxs = [3, 14]
    recognition_site = 'GAATTC'

    affected_codons = get_affected_codons_by_recognition_sites(dna_sequence, recognition_idxs, recognition_site)

    # Check that the correct affected codons are returned for each recognition site
    assert len(affected_codons) == 2  # There should be affected codons for both sites
    assert affected_codons[3] == [('GAA', [3, 4, 5]), ('TTC', [6, 7, 8])]  # Site at index 3
    assert affected_codons[14] == [('CTT', [12, 13, 14]), ('AAG', [15, 16, 17]),
                                   ('GAT', [18, 19, 20])]  # Site at index 14
