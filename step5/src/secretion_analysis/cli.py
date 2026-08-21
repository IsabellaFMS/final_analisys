from __future__ import annotations

import argparse
from pathlib import Path

WORKFLOWS = ("full", "signalp", "deeptmhmm", "comparative")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the secretion analysis step."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Base directory used to resolve input and output paths.",
    )
    parser.add_argument(
        "--workflow",
        choices=WORKFLOWS,
        default="full",
        help="Workflow to execute.",
    )
    parser.add_argument("--signalp-txt", type=Path, default=None)
    parser.add_argument("--signalp-gff", type=Path, default=None)
    parser.add_argument("--deeptmhmm-gff", type=Path, default=None)
    parser.add_argument("--conservation-map", type=Path, default=None)
    parser.add_argument("--out-signalp", type=Path, default=None)
    parser.add_argument("--out-deep", type=Path, default=None)
    parser.add_argument("--out-merge", type=Path, default=None)
    return parser
