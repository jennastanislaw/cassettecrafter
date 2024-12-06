import sys
import os
import pytest
# Tests
# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from mixed_base_rec_site_check_utils import (
    expand_dna_sequence,
    find_non_canonical_bases,
    check_non_canonical_in_rec_sites,
    check_recognition_sites_in_expanded_sequences,
    append_valid_sequences
)

class TestExpandDnaSequence:
    def test_single_base(self):
        """Test single valid IUPAC codes."""
        assert expand_dna_sequence("A") == ["A"]
        assert expand_dna_sequence("N") == ["A", "C", "G", "T"]
        assert expand_dna_sequence("R") == ["A", "G"]

    def test_multiple_bases(self):
        """Test sequences with multiple valid IUPAC codes."""
        assert expand_dna_sequence("AC") == ["AC"]
        assert sorted(expand_dna_sequence("AN")) == ["AA", "AC", "AG", "AT"]
        assert sorted(expand_dna_sequence("NR")) == ["AA", "AG", "CA", "CG", "GA", "GG", "TA", "TG"]

    def test_invalid_base(self):
        """Test sequences with invalid bases."""
        with pytest.raises(ValueError, match="Invalid base 'X' in DNA sequence."):
            expand_dna_sequence("AX")

    def test_empty_sequence(self):
        """Test an empty DNA sequence."""
        assert expand_dna_sequence("") == [""]

    def test_complex_sequence(self):
        """Test a complex DNA sequence with multiple IUPAC codes."""
        dna_sequence = "RYN"
        expected = [
            "ACA", "ACC", "ACG", "ACT",
            "ATA", "ATC", "ATG", "ATT",
            "GCA", "GCC", "GCG", "GCT",
            "GTA", "GTC", "GTG", "GTT",
        ]
        assert sorted(expand_dna_sequence(dna_sequence)) == sorted(expected)

class TestFindNonCanonicalBases:
    def test_no_non_canonical_bases(self):
        """Test a sequence with only canonical bases."""
        sequence = "ACGTACGT"
        expected = []
        assert find_non_canonical_bases(sequence) == expected

    def test_all_non_canonical_bases(self):
        """Test a sequence with only non-canonical bases."""
        sequence = "RYWSKMBDHVN"
        expected = [(0, 'R'), (1, 'Y'), (2, 'W'), (3, 'S'), (4, 'K'),
                    (5, 'M'), (6, 'B'), (7, 'D'), (8, 'H'), (9, 'V'), (10, 'N')]
        assert find_non_canonical_bases(sequence) == expected

    def test_mixed_sequence(self):
        """Test a sequence with a mix of canonical and non-canonical bases."""
        sequence = "ACGTNACGTRY"
        expected = [(4, 'N'), (9, 'R'), (10, 'Y')]
        assert find_non_canonical_bases(sequence) == expected

    def test_empty_sequence(self):
        """Test an empty sequence."""
        sequence = ""
        expected = []
        assert find_non_canonical_bases(sequence) == expected

    def test_lowercase_bases(self):
        """Test a sequence with lowercase bases."""
        sequence = "acgtnRY"
        expected = [(4, 'n'), (5, 'R'), (6, 'Y')]
        assert find_non_canonical_bases(sequence) == expected

    def test_special_characters(self):
        """Test a sequence with special characters."""
        sequence = "ACGT@!NRY"
        expected = [(4, '@'), (5, '!'), (6, 'N'), (7, 'R'), (8, 'Y')]
        assert find_non_canonical_bases(sequence) == expected



class TestCheckNonCanonicalInRecSites:

    def test_no_non_canonical_bases(self):
        """Test with no non-canonical bases."""
        indices = []
        fwd_sites = [5, 15]
        rev_sites = [25]
        enzyme_oh_length = 6
        assert not check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)

    def test_non_canonical_outside_sites(self):
        """Test with non-canonical bases outside recognition sites."""
        indices = [(3, 'R'), (20, 'Y')]
        fwd_sites = [5, 10]
        rev_sites = [25]
        enzyme_oh_length = 6
        assert not check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)

    def test_non_canonical_in_fwd_site(self):
        """Test with a non-canonical base overlapping a forward recognition site."""
        indices = [(6, 'R')]
        fwd_sites = [5, 15]
        rev_sites = [25]
        enzyme_oh_length = 6
        assert check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)

    def test_non_canonical_in_rev_site(self):
        """Test with a non-canonical base overlapping a reverse recognition site."""
        indices = [(27, 'Y')]
        fwd_sites = [5, 15]
        rev_sites = [25]
        enzyme_oh_length = 6
        assert check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)

    def test_multiple_recognition_sites(self):
        """Test with multiple recognition sites."""
        indices = [(16, 'R'), (30, 'Y')]
        fwd_sites = [5, 15]
        rev_sites = [25, 35]
        enzyme_oh_length = 6
        assert check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)

    def test_edge_case_overlap_end_site(self):
        """Test with a non-canonical base at the end of a recognition site."""
        indices = [(9, 'R')]
        fwd_sites = [5]
        rev_sites = [25]
        enzyme_oh_length = 5
        assert check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)

    def test_edge_case_no_overlap_end_site(self):
        """Test with a non-canonical base at the end of a recognition site."""
        indices = [(10, 'R')]
        fwd_sites = [5]
        rev_sites = [25]
        enzyme_oh_length = 5
        assert not check_non_canonical_in_rec_sites(indices, fwd_sites, rev_sites, enzyme_oh_length)


class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, oh_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.OH_length = oh_length

    def __repr__(self):
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, oh_length={self.oh_length})"


class TestCheckRecognitionSitesInExpandedSequences:
    @pytest.fixture
    def enzyme(self):
        """Create a test enzyme."""
        return Enzyme("TestEnzyme", "GAATTC", "CTTAAG", 0, 6)

    def test_no_non_canonical_overlap(self, enzyme):
        """Test when no non-canonical bases overlap with recognition sites."""
        enzyme = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        indices = [(5, 'N')]
        expanded_sequences = ["ATGCATGCAT", "GCTAGCTAAG"]

        result = check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme)

        # Assertions
        assert result == [False, False]

    def test_non_canonical_overlap(self, enzyme):
        """Test when a non-canonical base overlaps with a recognition site."""
        enzyme = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        indices = [(5, 'N')]
        expanded_sequences = ["ATGCGAATTC", "GCTACTTAAG"]

        result = check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme)

        # Assertions
        assert result == [True, True]

    def test_multiple_non_canonical_bases(self, enzyme):
        """Test with multiple non-canonical bases."""
        enzyme = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        indices = [(5, 'N'), (8, 'R')]
        expanded_sequences = ["ATGCGAATTC", "ATGCCGTAAG"]

        result = check_recognition_sites_in_expanded_sequences(indices, expanded_sequences, enzyme)

        # Assertions
        assert result == [True, False]

class TestAppendValidSequences:
    @pytest.fixture
    def setup_data(self):
        """Fixture to provide common test data."""
        name = "TestSequence"
        seq = "ATGCTN"
        expanded_sequences = ["ATGCTA", "ATGCTG", "ATGCTC"]
        return name, seq, expanded_sequences

    def test_with_recognition_sites(self, setup_data):
        """Test appending expanded sequences when recognition sites exist."""
        name, seq, expanded_sequences = setup_data
        rec_sites = [True, True, True]
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions
        assert len(new_rows) == 3
        for i, row in enumerate(new_rows):
            assert row['name'] == name
            assert row['DNA'] == seq
            assert row['valid_mixed_bases'] == expanded_sequences[i]

    def test_without_recognition_sites(self, setup_data):
        """Test appending the original sequence when no recognition sites exist."""
        name, seq, expanded_sequences = setup_data
        rec_sites = [False, False, False]
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions
        assert len(new_rows) == 1
        assert new_rows[0]['name'] == name
        assert new_rows[0]['DNA'] == seq
        assert new_rows[0]['valid_mixed_bases'] == seq

    def test_mixed_recognition_sites(self, setup_data):
        """Test mixed recognition site results."""
        name, seq, expanded_sequences = setup_data
        rec_sites = [True, False, True]
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions
        assert len(new_rows) == 3
        assert new_rows[0]['valid_mixed_bases'] == "ATGCTA"
        assert new_rows[1]['valid_mixed_bases'] == "ATGCTG"
        assert new_rows[2]['valid_mixed_bases'] == "ATGCTC"

    def test_empty_expanded_sequences(self):
        """Test with no expanded sequences provided."""
        name = "TestSequence"
        seq = "ATGCTN"
        expanded_sequences = []
        rec_sites = []
        new_rows = []

        append_valid_sequences(name, seq, expanded_sequences, rec_sites, new_rows)

        # Assertions
        assert len(new_rows) == 1
        assert new_rows[0]['name'] == name
        assert new_rows[0]['DNA'] == seq
        assert new_rows[0]['valid_mixed_bases'] == seq