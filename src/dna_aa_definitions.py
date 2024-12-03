# Define the standard genetic code for DNA
CODON_TABLE_DNA = {
    'F': ['TTT', 'TTC'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'],
    'M': ['ATG'],  # Start codon
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGU', 'AGC'], # missing last two in old code
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAC', 'TAT'],
    'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'],
    'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    '*': ['TAA', 'TAG', 'TGA'],  # Stop codons
}

# Invert the codon table for easier lookup
CODON_TO_AMINO_ACID_DNA = {codon: aa for aa, codons in CODON_TABLE_DNA.items() for codon in codons}

MIXED_BASES = {
    'R': ['A', 'G'],  # Purines
    'Y': ['C', 'T'],  # Pyrimidines
    'S': ['C', 'G'],  # Strong (G or C)
    'W': ['A', 'T'],  # Weak (A or T)
    'K': ['G', 'T'],  # Keto (G or T)
    'M': ['A', 'C'],  # Amino (A or C)
}

MIXED_BASES_COMBO_TO_BASE = {"".join(val): key for key, val in MIXED_BASES.items()}
