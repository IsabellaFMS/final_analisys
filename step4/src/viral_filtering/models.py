from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def protein_id(self) -> str:
        return self.header.replace("|", " ").split()[0]


@dataclass(frozen=True)
class ExtractionSummary:
    raw_rows: int
    valid_rows: int
    exclusion_ids: set[str] = field(default_factory=set)

    @property
    def unique_count(self) -> int:
        return len(self.exclusion_ids)


@dataclass(frozen=True)
class FilteringSummary:
    total_records: int
    excluded_records: int
    safe_records: list[FastaRecord] = field(default_factory=list)

    @property
    def safe_count(self) -> int:
        return len(self.safe_records)


@dataclass(frozen=True)
class DuplicateVerdict:
    protein_id: str
    kept_instead: str
    pathogen: str
    twin_excluded: bool
    source_index: int


@dataclass(frozen=True)
class ReincorporationSummary:
    verdicts: list[DuplicateVerdict] = field(default_factory=list)
    safe_records: list[FastaRecord] = field(default_factory=list)

    @property
    def reincorporated_count(self) -> int:
        return len(self.verdicts)

    @property
    def twin_excluded_count(self) -> int:
        return sum(1 for verdict in self.verdicts if verdict.twin_excluded)

    @property
    def safe_count(self) -> int:
        return len(self.safe_records)


@dataclass(frozen=True)
class StepResult:
    extraction: ExtractionSummary | None = None
    filtering: FilteringSummary | None = None
    reincorporation: ReincorporationSummary | None = None
