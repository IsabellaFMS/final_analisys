from __future__ import annotations

from pathlib import Path

from .models import DuplicateVerdict, FastaRecord


def _prepare_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_exclusion_ids(exclusion_ids: set[str], output_path: Path) -> None:
    _prepare_parent(output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        for value in sorted(exclusion_ids):
            handle.write(f"{value}\n")


def write_fasta_records(records: list[FastaRecord], output_path: Path) -> None:
    _prepare_parent(output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n{record.sequence}\n")


def write_flagged_verdicts(verdicts: list[DuplicateVerdict], output_path: Path) -> int:
    # The reincorporated accessions whose kept twin was excluded for host
    # homology share every 12-mer with it, so they carry the same homology and
    # the deposit has to name them.
    flagged = sorted((verdict for verdict in verdicts if verdict.twin_excluded), key=lambda verdict: verdict.source_index)
    _prepare_parent(output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("protein_id\tkept_instead\tpathogen\tsource_index\n")
        for verdict in flagged:
            handle.write(f"{verdict.protein_id}\t{verdict.kept_instead}\t{verdict.pathogen}\t{verdict.source_index}\n")
    return len(flagged)
