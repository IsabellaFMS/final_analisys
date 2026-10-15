from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class StepPaths:
    root: Path
    outputs: Path
    tables: Path
    reports: Path
    bcell: Path


@dataclass(frozen=True)
class Orthogroup:
    identifier: str
    members: tuple[str, ...]
    species: frozenset[str]
    consensus_length: float

    @property
    def is_conserved(self) -> bool:
        return len(self.species) == 4


@dataclass(frozen=True)
class Candidate:
    """A protein with everything the five rules need to decide on it."""

    accession: str
    description: str
    species: str
    orthogroup: str
    length: int
    consensus_length: float
    accessibility_group: str
    epitopes_high: int
    epitopes_intermediate: int
    mhc_class_i_alleles: int
    mhc_class_ii_alleles: int
    partial_description: bool

    @property
    def length_ratio(self) -> float:
        if not self.consensus_length:
            return 0.0
        return self.length / self.consensus_length

    @property
    def allele_coverage(self) -> int:
        return self.mhc_class_i_alleles + self.mhc_class_ii_alleles


@dataclass(frozen=True)
class Selection:
    candidate: Candidate
    stringency: float
    epitope_count: int
    replaces: str | None = None


@dataclass(frozen=True)
class TraceEntry:
    accession: str
    orthogroup: str
    verdict: str
    detail: str


@dataclass(frozen=True)
class SelectionResult:
    selected: list[Selection] = field(default_factory=list)
    trace: list[TraceEntry] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)
    parameters: dict[str, object] = field(default_factory=dict)
