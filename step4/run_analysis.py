from __future__ import annotations

import argparse
import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from viral_filtering import format_summary, resolve_command_paths, run_extraction, run_filtering, run_full_step


def _add_reincorporation_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--flagged-output", help="Flagged reincorporated accessions. Defaults to outputs/reports/flagged_reintegrated_accessions.txt.")
    parser.add_argument("--duplicate-report", help="Step 3 removed-records CSV. Defaults to ../step3/outputs/reports/removed_viral_proteins.csv.")
    parser.add_argument("--duplicate-fasta", help="Step 3 removed-records FASTA. Defaults to ../step3/outputs/fasta/removed_viral_proteins.fasta.")


def _add_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--tsv-input", help="TSV hit matrix. Defaults to sources/viral_hits.tsv.")
    parser.add_argument("--fasta-input", help="Viral FASTA input. Defaults to sources/clean_viral_proteins.fasta.")
    parser.add_argument("--exclusion-output", help="Exclusion list output. Defaults to outputs/reports/exclusion_ids.txt.")
    parser.add_argument("--fasta-output", help="Filtered FASTA output. Defaults to outputs/fasta/safe_sequences.fasta.")
    _add_reincorporation_arguments(parser)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the viral hit extraction and FASTA filtering step.")
    parser.add_argument(
        "--project-root",
        default=Path(__file__).resolve().parent,
        help="Base directory used to resolve input and output paths.",
    )
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser("run", help="Extract exclusion IDs and filter the FASTA in one workflow.")
    _add_run_arguments(run_parser)

    extract_parser = subparsers.add_parser("extract", help="Generate the exclusion list from the TSV hit matrix.")
    extract_parser.add_argument("--tsv-input", help="TSV hit matrix. Defaults to sources/viral_hits.tsv.")
    extract_parser.add_argument("--exclusion-output", help="Exclusion list output. Defaults to outputs/reports/exclusion_ids.txt.")

    filter_parser = subparsers.add_parser("filter", help="Filter a FASTA using an exclusion list.")
    filter_parser.add_argument("--fasta-input", help="Viral FASTA input. Defaults to sources/clean_viral_proteins.fasta.")
    filter_parser.add_argument("--exclusion-input", help="Exclusion list input. Defaults to outputs/reports/exclusion_ids.txt.")
    filter_parser.add_argument("--fasta-output", help="Filtered FASTA output. Defaults to outputs/fasta/safe_sequences.fasta.")
    _add_reincorporation_arguments(filter_parser)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    command = args.command or "run"
    try:
        paths = resolve_command_paths(command=command, project_root=Path(args.project_root), args=args)
        if command == "extract":
            result = run_extraction(paths)
        elif command == "filter":
            result = run_filtering(paths)
        else:
            result = run_full_step(paths)
        print(format_summary(command, result, paths.flagged_output))
        return 0
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
