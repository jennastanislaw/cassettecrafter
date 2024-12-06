import sys
import os
import pytest

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
from enzyme_site_replacement import replace_enzyme_site


class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, oh_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = spacer_length
        self.oh_length = oh_length

    def __repr__(self):
        return f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, oh_length={self.oh_length})"


class TestReplaceEnzymeSite:
    @pytest.fixture
    def sample_enzyme_data(self):
        """Create a list of demo enzyme objects."""
        EcoRI = Enzyme("EcoRI", "GAATTC", "CTTAAG", 0, 4)
        BamHI = Enzyme("BamHI", "GGATCC", "CCTAGG", 0, 4)
        FakeEnzyme = Enzyme("FakeEnzyme", "TGG", "TGG", 0, 4)

        return [EcoRI, BamHI, FakeEnzyme]

    def test_successful_replacement(self, sample_enzyme_data):
        """Test that a codon is successfully replaced for a recognition site."""
        # Extract EcoRI enzyme from sample data
        EcoRI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "EcoRI")
        DNA = "AGCTGAATTCAGCTCTTAAGCTAGATG"  # DNA sequence with EcoRI recognition site
        new_DNA = replace_enzyme_site(EcoRI, DNA)

        # Check that the new DNA does not contain the original recognition site
        assert "GAATTC" not in new_DNA
        assert "GAA" in new_DNA or "TTC" in new_DNA  # Check for synonymous codons

    def test_no_replacement_possible(self, sample_enzyme_data):
        """Test the scenario where no replacement is possible."""
        # Extract FakeEnzyme from sample data
        FakeEnzyme = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "FakeEnzyme")
        DNA = "TGGTGAGATCCAGCTCTTAAGCTAGATG"  # DNA sequence with BamHI recognition site
        key = "0"

        # The recognition site should still be present in the new DNA
        with pytest.raises(ValueError, match=f"No suitable replacement found for recognition site at index {key}."):
            replace_enzyme_site(FakeEnzyme, DNA)

    def test_multiple_sites_replacement(self, sample_enzyme_data):
        """Test that multiple recognition sites can be replaced."""
        # Extract EcoRI enzyme from sample data
        EcoRI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "EcoRI")
        DNA = "AGCTGAATTCGGATCCGCTCTTAAGCTAGATG"  # Contains both EcoRI and BamHI sites
        new_DNA = replace_enzyme_site(EcoRI, DNA)

        # Ensure EcoRI site is replaced
        assert "GAATTC" not in new_DNA

        # Extract BamHI enzyme from sample data
        BamHI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "BamHI")
        new_DNA = replace_enzyme_site(BamHI, new_DNA)

        # Ensure BamHI site is replaced
        assert "GGATCC" not in new_DNA

    def test_empty_dna_sequence(self, sample_enzyme_data):
        """Test the behavior with an empty DNA sequence."""
        # Extract EcoRI enzyme from sample data
        EcoRI = next(enzyme for enzyme in sample_enzyme_data if enzyme.name == "EcoRI")
        DNA = ""
        new_DNA = replace_enzyme_site(EcoRI, DNA)

        # The output should still be an empty string
        assert new_DNA == ""
