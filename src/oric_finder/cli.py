"""Command-line interface for Ori-C Finder."""

import argparse
import sys

from .io import read_fasta, write_results
from .gc_skew import (
    calculate_gc_skew,
    calculate_cumulative_gc_skew,
    find_skew_regions,
)
from .motifs import find_dnaa_motifs
from .candidates import identify_candidate_regions, merge_regions
from .scoring import rank_candidates


def build_parser():
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="oric-finder",
        description="Predict bacterial replication origin (OriC) regions."
    )

    parser.add_argument(
        "genome",
        help="Input genome FASTA file."
    )

    parser.add_argument(
        "-o",
        "--output",
        default="oric_results.txt",
        help="Output file for predicted OriC candidates "
             "(default: oric_results.txt)."
    )

    parser.add_argument(
        "--window",
        type=int,
        default=1000,
        help="Window size for GC-skew calculation "
             "(default: 1000)."
    )

    parser.add_argument(
        "--step",
        type=int,
        default=100,
        help="Step size for GC-skew calculation "
             "(default: 100)."
    )

    return parser


def run_pipeline(sequence, window, step):
    """Run the complete OriC prediction pipeline."""

    # 1. Calculate GC skew
    skew = calculate_gc_skew(
        sequence,
        window=window,
        step=step
    )

    # 2. Calculate cumulative GC skew
    cumulative_skew = calculate_cumulative_gc_skew(sequence)

    # 3. Identify GC-skew candidate regions
    skew_regions = find_skew_regions(
        skew,
        cumulative_skew
    )

    # 4. Find DnaA-box-like motifs
    motifs = find_dnaa_motifs(sequence)

    # 5. Combine GC-skew and motif evidence
    candidates = identify_candidate_regions(
        skew_regions,
        motifs,
        len(sequence)
    )

    # 6. Merge overlapping/nearby candidates
    candidates = merge_regions(candidates)

    # 7. Score and rank candidates
    ranked_candidates = rank_candidates(candidates)

    return ranked_candidates


def main():
    """Main entry point for the oric-finder command."""

    parser = build_parser()
    args = parser.parse_args()

    try:
        # Read genome
        sequence = read_fasta(args.genome)

        # Run prediction
        candidates = run_pipeline(
            sequence,
            window=args.window,
            step=args.step
        )

        # Write results
        write_results(
            candidates,
            args.output
        )

        print(f"OriC prediction completed.")
        print(f"Genome: {args.genome}")
        print(f"Candidates: {len(candidates)}")
        print(f"Results: {args.output}")

    except FileNotFoundError:
        print(
            f"Error: input file not found: {args.genome}",
            file=sys.stderr
        )
        sys.exit(1)

    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)

    except Exception as error:
        print(
            f"Error: {error}",
            file=sys.stderr
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
