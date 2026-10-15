from __future__ import annotations

from pathlib import Path

from .models import FastaRecord


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
