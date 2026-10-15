from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


PATHOGEN_PATTERN = re.compile(r"\[(.*?)\]")

UNIPROT_DATABASES = frozenset({"sp", "tr"})

DUPLICATE_REASON = "Duplicate"
SHORT_REASON = "Short"


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def protein_id(self) -> str:
        """The accession, taken from the first field of the header."""
        first_field = self.header.split()[0]
        parts = first_field.split("|")
        if len(parts) >= 3 and parts[0] in UNIPROT_DATABASES:
            return parts[1]
        return self.header.replace("|", " ").split()[0]

    @property
    def pathogen(self) -> str:
        match = PATHOGEN_PATTERN.search(self.header)
        return match.group(1) if match else ""


@dataclass(frozen=True)
class RemovedRecord:
    record: FastaRecord
    reason: str
    kept_instead: str = ""
    source_index: int = -1


@dataclass(frozen=True)
class ProteomeCleaningSummary:
    label: str
    source: Path
    total_read: int
    duplicates: int
    too_short: int
    kept: int
    repeated_ids: tuple[str, ...] = ()

    @property
    def removed(self) -> int:
        return self.duplicates + self.too_short


@dataclass(frozen=True)
class StepResult:
    kmer_size: int
    human: ProteomeCleaningSummary | None = None
    viral: ProteomeCleaningSummary | None = None

    @property
    def summaries(self) -> tuple[ProteomeCleaningSummary, ...]:
        return tuple(item for item in (self.human, self.viral) if item is not None)
