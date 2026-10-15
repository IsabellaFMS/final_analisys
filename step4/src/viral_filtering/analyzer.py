from __future__ import annotations

from collections.abc import Iterable

from .loaders import DUPLICATE_REASON
from .models import DuplicateVerdict, ExtractionSummary, FastaRecord, FilteringSummary, ReincorporationSummary


def extract_exclusion_ids(query_names: Iterable[str | None]) -> ExtractionSummary:
    # Collect the accessions of the viral proteins that matched the host.
    exclusion_ids: set[str] = set()
    raw_rows = 0
    valid_rows = 0
    for query_name in query_names:
        raw_rows += 1
        if query_name is None:
            continue
        cleaned = query_name.split("|", 1)[0].strip()
        if not cleaned:
            continue
        valid_rows += 1
        exclusion_ids.add(cleaned)
    return ExtractionSummary(raw_rows=raw_rows, valid_rows=valid_rows, exclusion_ids=exclusion_ids)


def filter_records_by_exclusion(records: Iterable[FastaRecord], exclusion_ids: set[str]) -> FilteringSummary:
    safe_records: list[FastaRecord] = []
    total_records = 0
    excluded_records = 0
    for record in records:
        total_records += 1
        if record.protein_id in exclusion_ids:
            excluded_records += 1
        else:
            safe_records.append(record)
    return FilteringSummary(total_records=total_records, excluded_records=excluded_records, safe_records=safe_records)


def _raw_positions(clean_count: int, removed_indices: set[int]) -> list[int]:
    # The cleaned FASTA is the raw one minus the removed records, in the same
    # order, so the n-th kept record sits at the n-th index step 3 did not drop.
    positions: list[int] = []
    index = 0
    while len(positions) < clean_count:
        if index not in removed_indices:
            positions.append(index)
        index += 1
    return positions


def reincorporate_duplicates(
    clean_records: list[FastaRecord],
    exclusion_ids: set[str],
    duplicate_rows: list[dict[str, str]],
    removed_records: list[FastaRecord],
) -> ReincorporationSummary:
    # Step 3 dropped sequences byte-identical to one it kept. They carry their
    # own accession, species and conservation class, so they go back into the
    # filtered proteome at the position they held in the raw file.
    removed_indices = {int(row["source_index"]) for row in duplicate_rows}
    expected_raw = len(clean_records) + len(duplicate_rows)
    if removed_indices and max(removed_indices) >= expected_raw:
        raise ValueError(
            "The duplicate report does not describe the FASTA that was filtered: "
            f"{len(clean_records)} records read, {len(duplicate_rows)} removed, "
            f"but the report points at index {max(removed_indices)}."
        )

    positions = _raw_positions(len(clean_records), removed_indices)
    placed = [
        (positions[index], record)
        for index, record in enumerate(clean_records)
        if record.protein_id not in exclusion_ids
    ]

    sequences = {record.protein_id: record for record in removed_records}
    verdicts: list[DuplicateVerdict] = []
    for row in duplicate_rows:
        if row.get("reason") != DUPLICATE_REASON:
            continue
        protein_id = row["protein_id"].strip()
        record = sequences.get(protein_id)
        if record is None:
            raise ValueError(f"The removed-records FASTA has no sequence for {protein_id}.")
        kept_instead = row.get("kept_instead", "").strip()
        verdicts.append(
            DuplicateVerdict(
                protein_id=protein_id,
                kept_instead=kept_instead,
                pathogen=row.get("pathogen", "").strip(),
                twin_excluded=kept_instead in exclusion_ids,
                source_index=int(row["source_index"]),
            )
        )
        placed.append((int(row["source_index"]), record))

    placed.sort(key=lambda item: item[0])
    return ReincorporationSummary(verdicts=verdicts, safe_records=[record for _, record in placed])
