from __future__ import annotations

from pathlib import Path

from .models import StepResult


def format_summary(command: str, result: StepResult, flagged_path: Path | None = None) -> str:
    lines = [
        "Execution summary",
        f"Workflow: {command}",
    ]
    if result.extraction is not None:
        lines.extend(
            [
                f"TSV rows scanned: {result.extraction.raw_rows}",
                f"Valid hit rows: {result.extraction.valid_rows}",
                f"Unique exclusions: {result.extraction.unique_count}",
            ]
        )
    if result.filtering is not None:
        lines.extend(
            [
                f"FASTA records scanned: {result.filtering.total_records}",
                f"Records excluded: {result.filtering.excluded_records}",
                f"Safe records after filtering: {result.filtering.safe_count}",
            ]
        )
    if result.reincorporation is not None:
        lines.extend(
            [
                f"Duplicates reincorporated: {result.reincorporation.reincorporated_count}",
                f"  of which identical to an excluded accession: {result.reincorporation.twin_excluded_count}",
                f"Safe records written: {result.reincorporation.safe_count}",
            ]
        )
        if flagged_path is not None:
            lines.append(f"Flagged accessions written to: {flagged_path}")
    return "\n".join(lines)
