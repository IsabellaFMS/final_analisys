from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from secretion_analysis.cli import build_parser
from secretion_analysis.config import resolve_step_paths
from secretion_analysis.reporting import format_execution_summary
from secretion_analysis.workflow import (
    run_comparative_workflow,
    run_deeptmhmm_workflow,
    run_full_workflow,
    run_signalp_workflow,
)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project_root = args.project_root.resolve() if args.project_root else PROJECT_ROOT
    paths = resolve_step_paths(
        project_root=project_root,
        workflow=args.workflow,
        signalp_txt=args.signalp_txt,
        signalp_gff=args.signalp_gff,
        deeptmhmm_gff=args.deeptmhmm_gff,
        conservation_map=args.conservation_map,
        signalp_output=args.out_signalp,
        deeptmhmm_output=args.out_deep,
        comparative_output=args.out_merge,
    )

    workflows = {
        "full": run_full_workflow,
        "signalp": run_signalp_workflow,
        "deeptmhmm": run_deeptmhmm_workflow,
        "comparative": run_comparative_workflow,
    }

    try:
        result = workflows[args.workflow](paths)
    except Exception as exc:  # pragma: no cover
        print(f"Execution failed: {exc}", file=sys.stderr)
        return 1

    print(format_execution_summary(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
