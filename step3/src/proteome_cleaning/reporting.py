from __future__ import annotations

from .models import ProteomeCleaningSummary, StepResult


def _format_proteome(summary: ProteomeCleaningSummary, kmer_size: int) -> list[str]:
    lines = [
        "",
        f"{summary.label}  ({summary.source.name})",
        f"  Records read:      {summary.total_read}",
        f"  Duplicates:        {summary.duplicates}",
        f"  Too short (<{kmer_size}):   {summary.too_short}",
        f"  Records removed:   {summary.removed}",
        f"  Records written:   {summary.kept}",
    ]
    if summary.repeated_ids:
        shown = ", ".join(summary.repeated_ids[:5])
        suffix = ", ..." if len(summary.repeated_ids) > 5 else ""
        lines.append(f"  WARNING repeated accessions: {len(summary.repeated_ids)} ({shown}{suffix})")
    return lines


def format_summary(command: str, result: StepResult) -> str:
    lines = [
        "Execution summary",
        f"Workflow: {command}",
        f"K-mer size: {result.kmer_size}",
    ]
    for summary in result.summaries:
        lines.extend(_format_proteome(summary, result.kmer_size))
    return "\n".join(lines)
