from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List


@dataclass
class ProteinRecord:
    accession: str
    residues: str = ""
    scores: List[float] = field(default_factory=list)


@dataclass
class EpitopeRecord:
    positions: List[int]
    residues: List[str]
    scores: List[float]

    @property
    def start_position(self) -> int:
        return self.positions[0] + 1 if self.positions else 0

    @property
    def end_position(self) -> int:
        return self.positions[-1] + 1 if self.positions else 0

    @property
    def length(self) -> int:
        return len(self.positions)

    @property
    def residue_sequence(self) -> str:
        return "".join(self.residues)


@dataclass(frozen=True)
class StepPaths:
    project_root: Path
    input_file: Path
    output_root: Path
    tables_dir: Path
    fasta_dir: Path
    conservation_map: Path


@dataclass(frozen=True)
class OutputArtifacts:
    epitope_table: Path
    fasta_files: List[Path]


@dataclass(frozen=True)
class ExecutionSummary:
    input_file: Path
    output_root: Path
    proteins_loaded: int
    accessions_with_epitopes: int
    epitope_count: int
    epitope_table: Path
    fasta_directory: Path
