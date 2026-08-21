from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

from .models import FastaRecord, RemovedRecord


BASE_COLUMNS = ("query_seqname", "protein_id", "reason", "sequence_length", "kept_instead", "source_index")
PATHOGEN_COLUMN = "pathogen"


def _prepare_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_fasta_records(records: Sequence[FastaRecord], output_path: Path) -> None:
    _prepare_parent(output_path)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for record in records:
            handle.write(f">{record.header}\n{record.sequence}\n")


def write_removed_report(
    removed: Sequence[RemovedRecord],
    output_path: Path,
    include_pathogen: bool,
) -> None:
    
    columns = list(BASE_COLUMNS)
    if include_pathogen:
        columns.append(PATHOGEN_COLUMN)
    _prepare_parent(output_path)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for entry in removed:
            row = [
                entry.record.header,
                entry.record.protein_id,
                entry.reason,
                len(entry.record.sequence),
                entry.kept_instead,
                entry.source_index,
            ]
            if include_pathogen:
                row.append(entry.record.pathogen)
            writer.writerow(row)


def write_summary_file(text: str, output_path: Path) -> None:
    _prepare_parent(output_path)
    output_path.write_text(f"{text}\n", encoding="utf-8", newline="")
