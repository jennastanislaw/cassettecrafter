from dna_aa_definitions import CODON_TABLE_DNA, CODON_TO_AMINO_ACID_DNA, MIXED_BASES, MIXED_BASES_COMBO_TO_BASE

def make_mut_dict(editable_codon_seq, all_combinations, name, mutations_df):
    """From a list of mutation combinations, generate a dicgtionary containing 
        mutation names, and their corresponding mutated sequence. The original
        sequence is 'mutated' by replacing the codon at the specified position
        with that codon that codes for one of the allowed amino acid mutations 
        at that position.

    Args:
        editable_codon_seq (list): list of codons, which begins as the codons of the
            original input gene sequence. Mutatable positions are modified, and 
            alternative codons are inserted at allowed positions
        all_combinations (list): List of tuples, containing all possible combinations
            of codons
        name (str): name of the gene whose sequences are being generated
        mutations_df (pandas DataFrame): dataframe containing data about allowed
            mutations at different positions in the original sequence

    Returns:
        dict: dictionary where keys are the names of mutants, and values are the 
            DNA sequence corresponding to that mutant
    """
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

            #Special case for naming mutations encoded by mixed base codons
            last_base = combo[i][-1]
            if last_base not in "ACGT":
                codon1 = combo[i][:2] + MIXED_BASES[last_base][0]
                codon2 = combo[i][:2] + MIXED_BASES[last_base][1]
                mut_aa = CODON_TO_AMINO_ACID_DNA[codon1]+"/"+CODON_TO_AMINO_ACID_DNA[codon2]
            else: 
                mut_aa = CODON_TO_AMINO_ACID_DNA[combo[i]]
            if og_aa != mut_aa:
                # If yes, add mutation info to name of sequence
                # Format: [og AA][position][mut AA] (example: H3A, mutate His at pos 3 to Ala)
                mut_name+=f"_{og_aa}{str(pos)}{mut_aa}"

        # Turn the list into a string - this is a possible DNA sequence
        seq="".join(editable_codon_seq)
        mut_library[mut_name] = seq
    
    return mut_library

def split_to_codons(dna_seq):
    """Split a DNA sequence to a list of codons by dividing it into fragments
        of length 3

    Args:
        dna_seq (str): string containing DNA sequence

    Returns:
        list: list of strings, where each item in the list corresponds to one 
            codon (set of 3 base pairs) from the input sequence
    """
    if type(dna_seq) != str:
        raise TypeError("DNA sequence must be a string") 
    elif len(dna_seq) % 3 != 0:
        # check that length of sequence is divisible by 3
        raise ValueError("Sequence must have a length that is a multiple of 3.")

    codon_list = [] # output list

    # Iterate over the sequence in steps of 3 to get codons
    for i in range(0, len(dna_seq), 3):
        codon = str(dna_seq[i:i + 3])
        codon_list.append(codon)
    
    return codon_list
