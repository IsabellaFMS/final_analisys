from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

from .models import DUPLICATE_REASON, SHORT_REASON, FastaRecord, RemovedRecord


def clean_records(
    records: Iterable[FastaRecord],
    kmer_size: int,
) -> tuple[list[FastaRecord], list[RemovedRecord]]:
    # Sequences shorter than k are dropped first, so a short record is never
    # picked as the survivor of a group of identical sequences.
    kept: list[FastaRecord] = []
    removed: list[RemovedRecord] = []
    first_seen: dict[str, FastaRecord] = {}
    for source_index, record in enumerate(records):
        if len(record.sequence) < kmer_size:
            removed.append(RemovedRecord(record=record, reason=SHORT_REASON, source_index=source_index))
            continue
        original = first_seen.get(record.sequence)
        if original is None:
            first_seen[record.sequence] = record
            kept.append(record)
        else:
            removed.append(
                RemovedRecord(
                    record=record,
                    reason=DUPLICATE_REASON,
                    kept_instead=original.protein_id,
                    source_index=source_index,
                )
            )
    return kept, removed


def find_repeated_ids(records: Iterable[FastaRecord]) -> tuple[str, ...]:
    # Two records under the same accession would make the step 4 match
    # ambiguous, so they are reported even though nothing is removed.
    counts = Counter(record.protein_id for record in records)
    return tuple(sorted(protein_id for protein_id, total in counts.items() if total > 1))
