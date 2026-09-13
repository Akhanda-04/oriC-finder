from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_fasta(filepath):
    """
    Read a FASTA file and return a SeqRecord.

    Parameters
    ----------
    filepath : str or Path
        Path to the FASTA file.

    Returns
    -------
    SeqRecord
        The genome sequence as a Biopython SeqRecord.

    Raises
    ------
    FileNotFoundError
        If the FASTA file does not exist.
    ValueError
        If the FASTA file is empty or contains multiple sequences.
    """

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    records = list(SeqIO.parse(filepath, "fasta"))

    if not records:
        raise ValueError(f"No sequences found in FASTA file: {filepath}")

    if len(records) > 1:
        raise ValueError(
            "Expected a single genome sequence, "
            f"but found {len(records)} sequences."
        )

    return records[0]


def write_results(results, filepath):
    """
    Write Ori-C candidate results to a tab-separated file.

    Parameters
    ----------
    results : list of dict
        Candidate results produced by the scoring/ranking pipeline.

    filepath : str or Path
        Output file path.
    """

    filepath = Path(filepath)

    # Create parent directory if it does not exist
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if not results:
        filepath.write_text("No Ori-C candidates found.\n")
        return

    # Use the keys from the first result as column names
    fields = list(results[0].keys())

    with filepath.open("w", encoding="utf-8") as output:
        # Header
        output.write("\t".join(fields) + "\n")

        # Results
        for result in results:
            row = []

            for field in fields:
                value = result.get(field, "")
                row.append(str(value))

            output.write("\t".join(row) + "\n")