from __future__ import annotations

from pathlib import Path

from .analyzer import extract_exclusion_ids, filter_records_by_exclusion, reincorporate_duplicates
from .config import CommandPaths
from .loaders import load_duplicate_report, load_exclusion_ids, load_fasta_records, load_query_names
from .models import ReincorporationSummary, StepResult
from .writer import write_exclusion_ids, write_fasta_records


def _filter_and_write(paths: CommandPaths, records, exclusion_ids: set[str]):
    filtering = filter_records_by_exclusion(records=records, exclusion_ids=exclusion_ids)
    reincorporation = _reincorporate(paths, records, exclusion_ids)
    write_fasta_records(reincorporation.safe_records, paths.fasta_output)
    return filtering, reincorporation


def _require(path: Path | None, label: str) -> Path:
    # Required, not optional: a missing step 3 artifact would silently cost the
    # proteome its 45 deduplicated accessions, and the step must say so instead.
    if path is None or not path.is_file():
        raise FileNotFoundError(
            f"{label} not found: {path}. It is written by step 3, so run step 3 "
            "first, or point at the file with the matching command-line option."
        )
    return path


def _reincorporate(paths: CommandPaths, records, exclusion_ids: set[str]) -> ReincorporationSummary:
    reincorporation = reincorporate_duplicates(
        clean_records=records,
        exclusion_ids=exclusion_ids,
        duplicate_rows=load_duplicate_report(_require(paths.duplicate_report, "Step 3 removal report")),
        removed_records=load_fasta_records(_require(paths.duplicate_fasta, "Step 3 removed-records FASTA")),
    )
    return reincorporation


def run_extraction(paths: CommandPaths) -> StepResult:
    extraction = extract_exclusion_ids(load_query_names(paths.tsv_input))
    write_exclusion_ids(extraction.exclusion_ids, paths.exclusion_output)
    return StepResult(extraction=extraction)


def run_filtering(paths: CommandPaths) -> StepResult:
    records = load_fasta_records(paths.fasta_input)
    filtering, reincorporation = _filter_and_write(paths, records, load_exclusion_ids(paths.exclusion_input))
    return StepResult(filtering=filtering, reincorporation=reincorporation)


def run_full_step(paths: CommandPaths) -> StepResult:
    extraction = extract_exclusion_ids(load_query_names(paths.tsv_input))
    write_exclusion_ids(extraction.exclusion_ids, paths.exclusion_output)
    records = load_fasta_records(paths.fasta_input)
    filtering, reincorporation = _filter_and_write(paths, records, extraction.exclusion_ids)
    return StepResult(extraction=extraction, filtering=filtering, reincorporation=reincorporation)
