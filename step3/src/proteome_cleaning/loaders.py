from __future__ import annotations

from pathlib import Path

from .models import FastaRecord


def load_fasta_records(fasta_path: Path) -> list[FastaRecord]:
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA input not found: {fasta_path}")
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_parts: list[str] = []
    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header=header, sequence="".join(sequence_parts)))
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(sequence_parts)))
    if not records:
        raise ValueError(f"FASTA input has no records: {fasta_path}")
    return records
