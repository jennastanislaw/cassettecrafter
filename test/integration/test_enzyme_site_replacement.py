import sys
import os
import pytest

# Add the src directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))
from enzyme_site_replacement import replace_enzyme_site

class TestReplaceEnzymeSite:
    @pytest.fixture
    def sample_enzyme_data(self, tmpdir):
        """Create a temporary CSV file with sample enzyme data."""
        csv_content = """Enzyme,Fwd_recognition_site,Rev_recognition_site,spacer_length,OH_length
EcoRI,GAATTC,CTTAAG,0,4
BamHI,GGATCC,CCTAGG,0,4
FakeEnzyme,TGG,TGG,0,4
"""
        csv_file = tmpdir.join("enzymes.csv")
        with open(csv_file, 'w') as f:
            f.write(csv_content)
        return str(csv_file)

    def test_successful_replacement(self, sample_enzyme_data):
        """Test that a codon is successfully replaced for a recognition site."""
        DNA = "AGCTGAATTCAGCTCTTAAGCTAGATG"  # DNA sequence with EcoRI recognition site
        enzyme_name = "EcoRI"
        new_DNA = replace_enzyme_site(sample_enzyme_data, enzyme_name, DNA)

        # Check that the new DNA does not contain the original recognition site
        assert "GAATTC" not in new_DNA
        assert "GAA" in new_DNA or "TTC" in new_DNA  # Check for synonymous codons

    def test_no_replacement_possible(self, sample_enzyme_data):
        """Test the scenario where no replacement is possible."""
        DNA = "TGGTGAGATCCAGCTCTTAAGCTAGATG"  # DNA sequence with BamHI recognition site
        enzyme_name = "FakeEnzyme"
        key="0"

        # The recognition site should still be present in the new DNA
        with pytest.raises (ValueError, match=f"No suitable replacement found for recognition site at index {key}."):
           replace_enzyme_site(sample_enzyme_data, enzyme_name, DNA)


    def test_multiple_sites_replacement(self, sample_enzyme_data):
        """Test that multiple recognition sites can be replaced."""
        DNA = "AGCTGAATTCGGATCCGCTCTTAAGCTAGATG"  # Contains both EcoRI and BamHI sites
        enzyme_name = "EcoRI"
        new_DNA = replace_enzyme_site(sample_enzyme_data, enzyme_name, DNA)

        # Ensure EcoRI site is replaced
        assert "GAATTC" not in new_DNA
        # Now test for BamHI
        enzyme_name = "BamHI"
        new_DNA = replace_enzyme_site(sample_enzyme_data, enzyme_name, new_DNA)

        # Ensure BamHI site is replaced
        assert "GGATCC" not in new_DNA

    def test_invalid_enzyme(self, sample_enzyme_data):
        """Test that an invalid enzyme name raises an appropriate error."""
        DNA = "AGCTGAATTCAGCTCTTAAGCTAGATG"
        enzyme_name = "InvalidEnzyme"

        with pytest.raises(KeyError):
            replace_enzyme_site(sample_enzyme_data, enzyme_name, DNA)

    def test_empty_dna_sequence(self, sample_enzyme_data):
        """Test the behavior with an empty DNA sequence."""
        DNA = ""
        enzyme_name = "EcoRI"
        new_DNA = replace_enzyme_site(sample_enzyme_data, enzyme_name, DNA)

        # The output should still be an empty string
        assert new_DNA == ""

