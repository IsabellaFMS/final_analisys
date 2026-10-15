#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from src.bcell_analysis import (
    ExtractionConfig,
    build_step_paths,
    format_execution_summary,
    run_bcell_analysis,
)

DEFAULT_CUTOFF = 0.1512
DEFAULT_INPUT_FILE = Path("sources/0.1512_raw_output.csv")


def build_parser() -> argparse.ArgumentParser:
    project_root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Run the B-cell epitope extraction step."
    )
    parser.add_argument(
        "--input-file",
        "-i",
        type=Path,
        default=DEFAULT_INPUT_FILE,
        help=f"Input CSV path, relative to the project root by default. Default: {DEFAULT_INPUT_FILE}.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path("outputs"),
        help="Output directory, relative to the project root by default.",
    )
    parser.add_argument(
        "--cutoff",
        "-c",
        type=float,
        default=DEFAULT_CUTOFF,
        help=f"Minimum residue score accepted as part of an epitope. Default: {DEFAULT_CUTOFF}.",
    )
    parser.add_argument(
        "--length-min",
        "-mi",
        dest="min_length",
        type=int,
        default=20,
        help="Minimum epitope length. Default: 20.",
    )
    parser.add_argument(
        "--length-max",
        "-ma",
        dest="max_length",
        type=int,
        default=30,
        help="Maximum epitope length. Default: 30.",
    )
    parser.add_argument(
        "--separator",
        "-s",
        default=",",
        help="CSV separator. Default: ','.",
    )
    parser.add_argument(
        "--classification-file",
        type=Path,
        default=None,
        help=(
            "Step 1 conservation report. Defaults to ../step1/outputs/ and falls "
            "back to the copy under sources/."
        ),
    )
    parser.epilog = (
        "This runner parses CLI arguments, resolves project-relative paths, "
        "calls reusable package APIs, and prints a concise summary. "
        f"Project root: {project_root}"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project_root = Path(__file__).resolve().parent

    config = ExtractionConfig(
        cutoff=args.cutoff,
        min_length=args.min_length,
        max_length=args.max_length,
        separator=args.separator,
    )
    paths = build_step_paths(
        project_root=project_root,
        input_file=args.input_file,
        output_dir=args.output_dir,
        classification_file=args.classification_file,
    )

    if not paths.input_file.exists():
        parser.error(
            f"Input file not found: {paths.input_file}. "
            f"Place the default file at '{DEFAULT_INPUT_FILE}' or pass --input-file."
        )

    if not paths.conservation_map.exists():
        parser.error(
            f"Conservation report not found: {paths.conservation_map}. "
            "It is written by step 1, so run step 1 first, or pass --classification-file."
        )

    summary = run_bcell_analysis(paths=paths, config=config)
    print(format_execution_summary(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
