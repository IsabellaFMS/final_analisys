from __future__ import annotations

import argparse
import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from proteome_cleaning import (
    format_summary,
    resolve_command_paths,
    run_full_step,
    run_human_cleaning,
    run_viral_cleaning,
)


def _add_kmer_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=12,
        help="Minimum sequence length to keep. Defaults to 12, the k used downstream.",
    )


def _add_human_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--human-input", help="Raw human FASTA. Defaults to sources/human.fasta.")
    parser.add_argument("--human-output", help="Cleaned human FASTA. Defaults to outputs/fasta/clean_human_proteome.fasta.")
    parser.add_argument("--human-report", help="Removed human records. Defaults to outputs/reports/removed_human_proteins.csv.")


def _add_viral_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--viral-input", help="Raw viral FASTA. Defaults to sources/viral_proteins.fasta.")
    parser.add_argument("--viral-output", help="Cleaned viral FASTA. Defaults to outputs/fasta/clean_viral_proteins.fasta.")
    parser.add_argument("--viral-report", help="Removed viral records. Defaults to outputs/reports/removed_viral_proteins.csv.")
    parser.add_argument("--viral-removed-output", help="FASTA of the removed viral records. Defaults to outputs/fasta/removed_viral_proteins.fasta.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean the raw proteomes before the k-mer comparison.")
    parser.add_argument(
        "--project-root",
        default=Path(__file__).resolve().parent,
        help="Base directory used to resolve input and output paths.",
    )
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser("run", help="Clean both proteomes in one workflow.")
    _add_kmer_argument(run_parser)
    _add_human_arguments(run_parser)
    _add_viral_arguments(run_parser)
    run_parser.add_argument("--summary-output", help="Run summary. Defaults to outputs/reports/cleaning_summary.txt.")

    human_parser = subparsers.add_parser("clean-human", help="Clean the human proteome only.")
    _add_kmer_argument(human_parser)
    _add_human_arguments(human_parser)

    viral_parser = subparsers.add_parser("clean-viral", help="Clean the viral proteome only.")
    _add_kmer_argument(viral_parser)
    _add_viral_arguments(viral_parser)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    command = args.command or "run"
    try:
        paths = resolve_command_paths(command=command, project_root=Path(args.project_root), args=args)
        if command == "clean-human":
            result = run_human_cleaning(paths)
        elif command == "clean-viral":
            result = run_viral_cleaning(paths)
        else:
            result = run_full_step(paths)
        print(format_summary(command, result))
        return 0
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
