import copy
from Bio.Seq import Seq
from data.dna_aa_definitions import CODON_TABLE_DNA, CODON_TO_AMINO_ACID_DNA, MIXED_BASES, MIXED_BASES_COMBO_TO_BASE

### Helper functions for read_input ###
# Get name from te fasta or csv
def get_gene_name(lines,filetype):
    """Get gene name from the first column of csv, or the first line of a fasta

    Args:
        lines (list): list of strings corresponding to lines in an input file
        filetype (str): string containing the file type of an input file 
            (determined by its extension)

    Returns:
        str : name as a string
    """
    # If csv, take name from the first column
    if filetype == "csv":
        name=lines[0].split(",")[0]
    # If fasta, take names from the first row
    elif filetype == "fa" or filetype == "fasta":
        name=lines[0].replace(">","")
    else:
        # return None
        raise ValueError(f"Filetype {filetype} is not allowed. Valid input file types are csv or fa")

    #remove extra characters if needed
    return name.strip().replace('\ufeff', '') 

# Read in the seq from the file
def get_seq(lines,filetype):
    """Get DNA sequence from the second column of csv, or the remaining lines (lines after name) of a fasta

    Args:
        lines (list): list of strings corresponding to lines in an input file
        filetype (str): string containing the file type of an input file 
            (determined by its extension)

    Returns:
        BioPython Seq: Seq object containing information about the original DNA 
            sequence of the gene
    """

    if filetype == "csv":
        # If file is csv, seq should be in second column. Assuming no header
        seq=lines[0].split(",")[-1].strip()
    elif filetype == "fa":
        # If file is fasta, first line will contain information
        # all consecutive lines should contain sequence, but the 
        # sequence may be split over multiple lines
        seq=""
        for line in lines[1:]:
            seq+=line.strip()

    return seq

def biopython_seq_from_str(str_seq):
    """Accepts an amino acid or DNA sequence and converts it to a BipPython
        Seq object. If this is an amino acid sequence, it will first be 
        reverse transcribed in DNA (see convert_aa_to_dna for more details) 

    Args:
        str_seq (str): amino acid or DNA sequence 

    Returns:
        BioPython Seq object : the input sequence as BioPython Seq object
    """
    assert type(str_seq) == str, f"The sequence much be a string, not a(n) {type(str_seq).__name__}"
    
    #make it a biopython sequence for processing later
    str_seq = str_seq.upper()
    seq_obj = Seq(str_seq)

    bases={"A","C","G","T"}

    # If this is not a DNA sequence, convert it 
    if not set(str_seq).issubset(bases):
        dna_seq = convert_aa_to_dna(seq_obj._data)
        seq_obj = Seq(dna_seq)
    
    return seq_obj

def convert_aa_to_dna(prot_seq, codon_library=CODON_TABLE_DNA):
    """Converts and amino acid string to a DNA sequence composed of codons that
        correspond to the correct amino acid. The codons are selected by the 
        first corresponding item in the codon library.

    Args:
        prot_seq (str): string of amino acids, using the one-letter code
        codon_library (dictionary, optional): Pre-defined dictionary mapping
            amino acid keys to a list of possible codons. The first element in
            each list will be the codon used in the sequence that is returned.
            Defaults to CODON_TABLE_DNA.

    Raises:
        ValueError: if amino acid is invalid (not one of the 20 canonical AAs)

    Returns:
        str : the converted DNA sequence
    """
    assert type(prot_seq) == str, f"The sequence must be a string, not a(n) {type(prot_seq).__name__}"
    
    dna_seq = ""
    allowed_aas = list(codon_library.keys())
    for aa in prot_seq.upper():
        if aa in allowed_aas:
            dna_seq += codon_library[aa][0] # select the first codon for the amino acid
        else:
            raise ValueError(f"Invalid amino acid: {aa}")
    
    return dna_seq

# Take comma separated list and return split list of strings
def split_csl(csl):
    """Splits an input string that contains commas into a list of string, 
        where each item in the list correponds a value in the original string list

        For example, for the input 'A,B,C', the output would be ['A','B','C']

    Args:
        csl (str): "comma-separated-list"; string containing a list that is delimited
            by commas

    Returns:
        list: list of strings
    """
    if type(csl) != str:
        raise TypeError(f"Input must be a string, not a(n) {type(csl).__name__}")
    elif "," not in csl:
        raise ValueError("Input must be formatted as a comma-separated list")
   
    return list(val.strip() for val in csl.split(","))

def get_allowed_codon_list(mutations_aa, original_codons):
    """Convert amino acid information to codon information for each of the
        allowed mutations

    Args:
        mutations_aa (Pandas Series): columns from Pandas DataFrame which contains a list
            of allowed amino acid mutations for each row (position in the sequence)
        original_codons (Pandas Series): columns from Pandas DataFrame which contains 
            the original codon that the corresponding position

    Returns:
        list: list of allowed codons. The order of this list correlates with the
            order of positions given in the mutation csv and its dataframe. Note this
            includes the original codon sequence at that position, because it is
            also possible to have no mutation at a given position
    """
    output=list()

    # For each positions's possible mutations, convert the list of amino acids
    # to a list of codons
    for pos,allowed_aa_at_pos in enumerate(mutations_aa):
        allowed_codons_at_pos= list() 
        original_codon = original_codons[pos]
        original_aa = CODON_TO_AMINO_ACID_DNA[original_codon]
    
        # Convert each AA in list to the first codon in the translation table
        for aa_mut in allowed_aa_at_pos:
            if not aa_mut == original_aa: #Don't add mutation to list if it is the same as the original
                allowed_codons_at_pos.append(CODON_TABLE_DNA[aa_mut][0])

        # Append list of allowed mutation codons and the original codon to output
        allowed_codons_at_pos.append(original_codon)

        final_allowed_codon_list = add_mixed_bases_and_combine(allowed_codons_at_pos)
        output.append(final_allowed_codon_list)

    return output

def add_mixed_bases_and_combine(allowed_codons_at_pos):
    """Accepts an input that consists of codons mapping to allowed mutations and
        checks this list to see if there are any codons that can be combined
        into a codon with a mixed based. Codons can include a mixed base which 
        can encode for multiple bases, allowing for two codons that encode
        two different amino acids but only vary in sequence at the third base
        to be encoded by a single codon instead. The simplified list with
        mixed base codon is returned

    Args:
        allowed_codons_at_pos (list of str): List of codons, corresponging to
            the allowed codons at one position in the input sequence

    Returns:
        list: list of strings of codons of allowed mutations
    """
    modifiable_allowed_codon_list = copy.deepcopy(allowed_codons_at_pos)
    if type(allowed_codons_at_pos) != list:
        raise TypeError("Allowed codons must be a list of strings")
    if len(allowed_codons_at_pos) == 0: #Return empty list if not mutations
        return []
    # Return the single mutation if there is only one
    # Need to capture edge case here, otherwise will throw error later
    elif len(allowed_codons_at_pos) < 1: 
        return allowed_codons_at_pos
    for i, curr_codon in enumerate(allowed_codons_at_pos[:-1]):
        curr_codon_mod = check_leu_arg_ser(curr_codon)
        base_1_2 = curr_codon_mod[:2]

        # Compare beginning of current codon to those in the allowed codons list
        for j, match_codon in enumerate(allowed_codons_at_pos[i+1:]):
            match_codon_mod = check_leu_arg_ser(match_codon)
            match_first_two = match_codon_mod[:2]

            # If there is a match, find the mixed base that fits both codons
            if match_first_two == base_1_2: # add something for special cases: ARG (R) and LEU (L)
                mixed_base_codon = get_mixed_base_codon(curr_codon_mod, match_codon_mod)
                if mixed_base_codon[:2] != base_1_2:
                    mixed_base_codon = base_1_2 + mixed_base_codon[-1]

                # Remove the two codons and replace with the mixed-base codon
                modifiable_allowed_codon_list.remove(curr_codon)
                modifiable_allowed_codon_list.remove(match_codon)
                modifiable_allowed_codon_list.append(mixed_base_codon)

                break #exit the inner loop, because no other matches could exist
    return modifiable_allowed_codon_list 

def get_mixed_base_codon(codonA, codonB):
    """Find the appropriate codon containing a mixed base in the third position
        which satisfies both of the input codons. In other words, the input
        codons share the same two bases, but differ on the third base. This 
        function returns a codon with a third base that is a mixed base and 
        not a conventional amino acid. This specific mixed base should encode 
        for both of the third bases of the input codon.

    Args:
        codonA (str): three-letter codon of one of the amino acids
        codonB (str): three-letter codon of the other amino acid

    Returns:
        str: codon containing mixed based
    """
    for input_codon in [codonA, codonB]:
        if type(input_codon) != str:
            raise TypeError("Input codon must be a string")
        elif len(input_codon) != 3:
            raise ValueError("Input codon must be a string with 3 letters")
    
    # Sort the two bases alphabetically to generate the correct key of lookup
    base3_l = [codonA[-1],codonB[-1]]
    base3_l.sort()
    base3_joint_str = "".join(base3_l)

    # Access the pre-defined dictionary mapping bases to mixed bases, and
    # outputs value paired with the key generated above
    base3_mixed = MIXED_BASES_COMBO_TO_BASE[base3_joint_str]

    # Full codon with the third position as a mixed base
    base12 = codonA[:2]
    codon = base12 + base3_mixed

    return codon

def check_leu_arg_ser(codon):
    """Check the provided codon and sees if it codes for leucine,
        arginine, or serine. These 3 amino acids each have 6 codons that encode
        them; however, only 2/6 of these codons share the first 2 bases with 
        the codon for another amino acid. Therefore, to capture these edge cases
        and check if there is a possibility to use a mixed codon, this function 
        returns one of the codons which has a possible match in the first 2 
        bases to another amino acid's codons.

        (Example for explanation:
        Leucine codons: ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
        Phenylalanine codons: ['TTT', 'TTC']
        Leu must be encoded by TTA or TTG for it to share the first two
        bases with Phe. Therefore, we may need to modify the codon being used
        to check if there is a possibiity to use a mixed base which encodes two
        possible amino acids at the same position 
        )

    Args:
        codon (str): three letter codon 

    Returns:
        tuple of strings: (a codon whose first two bases are shared with another
                amino acid's codons, the first two bases of that codon as a str)
    """
    # Check input type
    if type(codon) != str:
        raise TypeError("Input codon must be a string")
    elif len(codon) != 3:
        raise ValueError("Input codon must be a string with 3 letters")
    
    # Assign the first possible codon with a possible match to another amino acid codon if the first two bases
    codon_base_1_2 = codon[:2]
    if codon_base_1_2 == "CT": # Leu
        new_codon = "TTA"
    elif codon_base_1_2 == "CG": # Arg
        new_codon = "AGA"
    elif codon_base_1_2 == "TC": # Ser
        new_codon = "AGC"
    else:
        new_codon = codon
    return new_codon
