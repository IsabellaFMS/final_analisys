from __future__ import annotations

from pathlib import Path

from .analyzer import clean_records, find_repeated_ids
from .config import CommandPaths
from .loaders import load_fasta_records
from .models import DUPLICATE_REASON, SHORT_REASON, ProteomeCleaningSummary, StepResult
from .reporting import format_summary
from .writer import write_fasta_records, write_removed_report, write_summary_file


def _clean_proteome(
    label: str,
    source: Path,
    fasta_output: Path,
    report_output: Path,
    kmer_size: int,
    include_pathogen: bool,
    removed_fasta_output: Path | None = None,
) -> ProteomeCleaningSummary:
    records = load_fasta_records(source)
    kept, removed = clean_records(records, kmer_size)
    write_fasta_records(kept, fasta_output)
    write_removed_report(removed, report_output, include_pathogen=include_pathogen)
    if removed_fasta_output is not None:
        # Step 4 needs the sequences back to decide what to reincorporate.
        write_fasta_records([entry.record for entry in removed], removed_fasta_output)
    return ProteomeCleaningSummary(
        label=label,
        source=source,
        total_read=len(records),
        duplicates=sum(1 for entry in removed if entry.reason == DUPLICATE_REASON),
        too_short=sum(1 for entry in removed if entry.reason == SHORT_REASON),
        kept=len(kept),
        repeated_ids=find_repeated_ids(kept),
    )


def run_human_cleaning(paths: CommandPaths) -> StepResult:
    summary = _clean_proteome(
        label="human",
        source=paths.human_input,
        fasta_output=paths.human_output,
        report_output=paths.human_report,
        kmer_size=paths.kmer_size,
        include_pathogen=False,
    )
    return StepResult(kmer_size=paths.kmer_size, human=summary)


def run_viral_cleaning(paths: CommandPaths) -> StepResult:
    summary = _clean_proteome(
        label="viral",
        source=paths.viral_input,
        fasta_output=paths.viral_output,
        report_output=paths.viral_report,
        kmer_size=paths.kmer_size,
        include_pathogen=True,
        removed_fasta_output=paths.viral_removed_output,
    )
    return StepResult(kmer_size=paths.kmer_size, viral=summary)


def run_full_step(paths: CommandPaths) -> StepResult:
    result = StepResult(
        kmer_size=paths.kmer_size,
        human=run_human_cleaning(paths).human,
        viral=run_viral_cleaning(paths).viral,
    )
    if paths.summary_output is not None:
        write_summary_file(format_summary("run", result), paths.summary_output)
    return result
