from __future__ import annotations

import csv
from pathlib import Path

from .models import FastaRecord


def load_query_names(tsv_path: Path) -> list[str | None]:
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV input not found: {tsv_path}")
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "query_seqname" not in reader.fieldnames:
            raise ValueError("TSV input must contain a 'query_seqname' column.")
        return [row.get("query_seqname") for row in reader]


def load_exclusion_ids(path: Path) -> set[str]:
    if not path.exists():
        raise FileNotFoundError(f"Exclusion list not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


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
    return records


DUPLICATE_REASON = "Duplicate"
DUPLICATE_REPORT_COLUMNS = ("protein_id", "reason", "kept_instead", "source_index")


def load_duplicate_report(path: Path) -> list[dict[str, str]]:
    # The step 3 cleaning report: what was removed, why, and from where.
    if not path.exists():
        raise FileNotFoundError(f"Duplicate report not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in DUPLICATE_REPORT_COLUMNS if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Duplicate report is missing the column(s): {', '.join(missing)}")
        return list(reader)
