#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from antigen_selection import default_paths, format_summary, run_selection


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Select the final antigens from the step 5, 6 and 7 tables.")
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Re-run the step 6 extraction even when the epitope tables are already present.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    paths = default_paths()
    for directory in (paths.tables, paths.reports, paths.bcell):
        directory.mkdir(parents=True, exist_ok=True)

    try:
        result, _ = run_selection(paths, refresh=args.refresh)
    except (FileNotFoundError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    summary = format_summary(result)
    (paths.reports / "selection_summary.txt").write_text(summary, encoding="utf-8")
    print(summary, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
