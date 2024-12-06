

def find_starts_of_consecutive_indices(conserved_indices, consecutive_count):
    """
    Finds the indices in the conserved index list that are the start of a specified number
    of consecutive conserved indices.

    Args:
        conserved_indices (list): A sorted list of conserved indices.
        consecutive_count (int): The number of consecutive indices to look for (default is 4).

    Returns:
        list: Indices in the conserved index list that mark the start of consecutive conserved indices.
    """
    starts = []

    for i in range(len(conserved_indices) - consecutive_count + 1):
        # Check if the next (consecutive_count - 1) indices are consecutive
        if all(conserved_indices[j] + 1 == conserved_indices[j + 1] for j in range(i, i + consecutive_count - 1)):
            starts.append(i)  # Store the index of the start of the sequence

    return starts


def adjust_oligo_lengths(min_oligo_len, max_oligo_len, recognition_site_len, spacer_bps, OH_len):
    """
    Adjusts oligo length constraints based on enzyme properties.

    Parameters:
        min_oligo_len (int): Minimum oligo length.
        max_oligo_len (int): Maximum oligo length.
        recognition_site_len (int): Length of the enzyme recognition site.
        spacer_bps (int): Number of base pairs in the spacer.
        OH_len (int): Overhang length of the enzyme.

    Returns:
        tuple: Adjusted minimum and maximum oligo lengths.

    Raises:
        ValueError: If the max oligo length is less than the length adjustment required for the enzyme.
    """
    # Calculate the length adjustment
    len_adj = 2 * (recognition_site_len + spacer_bps + OH_len)

    # Check if max oligo length is compatible with the enzyme
    if max_oligo_len <= len_adj:
        raise ValueError(
            f"The max oligo length ({max_oligo_len}) is incompatible with the selected enzyme. "
            f"It must be at least {len_adj+1} to accommodate the enzyme properties "
            f"(recognition site: {recognition_site_len}, spacer: {spacer_bps}, overhang: {OH_len})."
        )

    # Adjust oligo lengths
    true_min_oligo_len = max(1, min_oligo_len - len_adj)
    true_max_oligo_len = max_oligo_len - len_adj

    return true_min_oligo_len, true_max_oligo_len



def calculate_oligo_lengths(reference, split_indices, OH_len):
    """Calculate lengths of oligos based on split indices."""

    if any(index < 0 or index > len(reference) for index in split_indices):
        raise ValueError("One or more split indices are out of the range of the reference sequence length.")

    oligo_lengths = []

    if not split_indices:
        # No splits, return the full length of the reference
        oligo_lengths.append(len(reference))
    else:
        # Calculate oligo lengths normally
        oligo_lengths.append(split_indices[0] + OH_len)  # From start to first index
        oligo_lengths.extend(split_indices[i + 1] - split_indices[i] + OH_len for i in range(len(split_indices) - 1))
        oligo_lengths.append(len(reference) - split_indices[-1])  # From last index to end

    return oligo_lengths


def generate_cassettes(DNA_df, split_sites, enzyme):

    OH_len = enzyme.OH_length

    split_sites[-1]=split_sites[-1]-OH_len

    for i in range(len(split_sites)-1):

        start = split_sites[i]
        end = split_sites[i + 1]+OH_len

        DNA_df[f'Cassette {i+1}'] = DNA_df['rec_sites_removed'].apply(lambda x: x[start:end])

    return DNA_df


def is_valid_overhang(overhang):
    """Checks if an overhang meets the defined rules."""
    # Rule 2: Avoid palindromes
    if len(overhang) == 0:
        return True

    if len(overhang) > 1 and overhang == overhang[::-1]:
        return False

    # Rule 3: No overhangs with the same three nucleotides in a row
    if any(overhang[i] == overhang[i+1] == overhang[i+2] for i in range(len(overhang) - 2)):
        return False

    # Rule 5: Avoid overhangs with extreme GC content (0% or 100%)
    gc_content = (overhang.count('G') + overhang.count('C')) / len(overhang)
    if len(overhang) > 1 and  gc_content == 0 or gc_content == 1:
        return False

    return True


def hamming_distance(seq1, seq2):
    """Calculates the Hamming distance between two sequences."""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def find_constant_indices(reference, sequences):
    """
    Find indices of positions that are constant across all sequences
    relative to the reference.

    Args:
        reference (str): The reference sequence.
        sequences (list of str): List of sequences to compare against the reference.

    Returns:
        list: Indices of positions conserved across all sequences.
    """
    # Initialize conserved indices with all positions of the reference
    conserved_indices = set(range(len(reference)))

    for seq in sequences:
        # Find indices where reference and seq characters match
        matching_indices = {i for i, (ref_char, seq_char) in enumerate(zip(reference, seq)) if ref_char == seq_char}
        # Update conserved indices to keep only those consistent across all sequences
        conserved_indices &= matching_indices

    return sorted(conserved_indices)
