from Bio import Align

def align_to_reference(reference, sequence):
    """
    Align a sequence to the reference sequence using pairwise alignment.

    Args:
        reference (str): The reference sequence.
        sequence (str): The sequence to align to the reference.

    Returns:
        tuple: Aligned reference and sequence.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(reference, sequence)
    # Access the first alignment's target and query (aligned reference and sequence)
    aligned_ref = alignments[0].target
    aligned_seq = alignments[0].query
    return aligned_ref, aligned_seq


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
    ref_aligned = reference

    for seq in sequences:
        # Align current sequence to the reference
        aligned_ref, aligned_seq = align_to_reference(ref_aligned, seq)

        # Update conserved indices by checking aligned positions
        new_conserved_indices = set()
        for i, (ref_char, seq_char) in enumerate(zip(aligned_ref, aligned_seq)):
            # Only keep positions where the characters are identical and not gaps
            if ref_char == seq_char and ref_char != "-":
                new_conserved_indices.add(i)

        # Update conserved indices
        conserved_indices &= new_conserved_indices

        # Update the reference with gaps to keep alignment consistent
        ref_aligned = aligned_ref

    return sorted(conserved_indices)


def find_starts_of_consecutive_indices(conserved_indices, consecutive_count=4):
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


