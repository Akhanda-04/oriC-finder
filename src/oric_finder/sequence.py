# sequence.py

def validate_sequence(sequence):
    """
    Validate whether a DNA sequence contains only A, T, G, and C.

    Parameters
    ----------
    sequence : str
        DNA sequence to validate.

    Returns
    -------
    bool
        True if the sequence is valid, otherwise False.
    """
    sequence = sequence.upper()

    valid_bases = {"A", "T", "G", "C"}

    return all(base in valid_bases for base in sequence)


def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    Returns
    -------
    str
        Reverse-complemented DNA sequence.
    """
    sequence = sequence.upper()

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    if not validate_sequence(sequence):
        raise ValueError("Invalid DNA sequence. Use only A, T, G, and C.")

    return "".join(complement[base] for base in reversed(sequence))


def gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    Returns
    -------
    float
        GC content as a percentage.
    """
    sequence = sequence.upper()

    if not validate_sequence(sequence):
        raise ValueError("Invalid DNA sequence. Use only A, T, G, and C.")

    if len(sequence) == 0:
        raise ValueError("Sequence cannot be empty.")

    gc_count = sequence.count("G") + sequence.count("C")

    return (gc_count / len(sequence)) * 100
