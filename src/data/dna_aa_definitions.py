# Define the standard genetic code for DNA
# From python_codon_tables https://pypi.org/project/python-codon-tables/ 
# Optimized for E. coli
CODON_TABLE_DNA={
    '*': ['TAA', 'TAG', 'TGA'], # Stop codons
    'A': ['GCA', 'GCC', 'GCG', 'GCT'], 
    'C': ['TGC', 'TGT'], 
    'D': ['GAC', 'GAT'], 
    'E': ['GAA', 'GAG'], 
    'F': ['TTC', 'TTT'], 
    'G': ['GGA', 'GGC', 'GGG', 'GGT'], 
    'H': ['CAC', 'CAT'], 
    'I': ['ATA', 'ATC', 'ATT'], 
    'K': ['AAA', 'AAG'], 
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 
    'M': ['ATG'],  # Start codon
    'N': ['AAC', 'AAT'], 
    'P': ['CCA', 'CCC', 'CCG', 'CCT'], 
    'Q': ['CAA', 'CAG'], 
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 
    'T': ['ACA', 'ACC', 'ACG', 'ACT'], 
    'V': ['GTA', 'GTC', 'GTG', 'GTT'], 
    'W': ['TGG'], 
    'Y': ['TAC', 'TAT']
}

# Invert the codon table for easier lookup
CODON_TO_AMINO_ACID_DNA = {codon: aa for aa, codons in CODON_TABLE_DNA.items() for codon in codons}

# Mixed bases that can be used to encode for multiple bases with one 
# See https://www.idtdna.com/pages/products/custom-dna-rna/mixed-bases
MIXED_BASES = {
    'R': ['A', 'G'],  # Purines
    'Y': ['C', 'T'],  # Pyrimidines
    'S': ['C', 'G'],  # Strong (G or C)
    'W': ['A', 'T'],  # Weak (A or T)
    'K': ['G', 'T'],  # Keto (G or T)
    'M': ['A', 'C'],  # Amino (A or C)
}

MIXED_BASES_COMBO_TO_BASE = {"".join(val): key for key, val in MIXED_BASES.items()}
