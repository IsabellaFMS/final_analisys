from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from epitope_step import StepPaths, format_summary, run_workflows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the epitope analysis step.",
    )
    parser.add_argument(
        "--workflow",
        choices=("all", "mhci", "mhcii"),
        default="all",
        help="Workflow to execute. 'all' runs MHCI and MHCII in sequence.",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT_DIR,
        help="Project root containing sources/ and outputs/.",
    )
    parser.add_argument(
        "--mhci-cutoff",
        type=float,
        default=0.75,
        help="Rank cutoff used by the MHCI workflow.",
    )
    parser.add_argument(
        "--mhcii-cutoffs",
        type=float,
        nargs="+",
        default=(5.0, 7.5, 10.0),
        help="Rank cutoffs used by the MHCII workflow before consolidation.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    paths = StepPaths.from_root(args.root)
    workflows = ("mhci", "mhcii") if args.workflow == "all" else (args.workflow,)
    results = run_workflows(
        paths=paths,
        workflows=workflows,
        mhci_cutoff=args.mhci_cutoff,
        mhcii_cutoffs=tuple(args.mhcii_cutoffs),
    )
    print(format_summary(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
