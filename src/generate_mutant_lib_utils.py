from dna_aa_definitions import CODON_TO_AMINO_ACID_DNA

def make_mut_dict(editable_codon_seq, all_combinations, name, mutations_df):
    mut_library = dict()
    positions=mutations_df.index.tolist() # positions where mut can be make

    # Loop over each combination
    for combo in all_combinations:
        mut_name=name #start building up name for mutant

        for i,pos in enumerate(positions):
            # Make "mutation" by inserting codon at allowable position 
            editable_codon_seq[pos-1] = combo[i]

            # Check if this combo contains a mutation from the original amino acid
            og_aa = mutations_df.loc[pos,"original"]
            mut_aa = CODON_TO_AMINO_ACID_DNA[combo[i]]
            if og_aa != mut_aa:
                # If yes, add mutation info to name of sequence
                # Format: [og AA][position][mut AA]
                mut_name+=f"_{og_aa}{str(pos)}{mut_aa}"

        seq="".join(editable_codon_seq)
        mut_library[mut_name] = seq
    
    return mut_library

def split_to_codons(dna_seq):
    codon_list = [] # output list

    # check that length of sequence is divisible by 3
    assert len(dna_seq) % 3 == 0

    # Iterate over the sequence in steps of 3 to get codons
    for i in range(0, len(dna_seq), 3):
        codon = str(dna_seq[i:i + 3])
        codon_list.append(codon)
    
    return codon_list
