from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class StepPaths:
    project_root: Path
    workflow: str
    sources_dir: Path
    outputs_dir: Path
    docs_dir: Path
    signalp_txt: Path
    signalp_gff: Path
    deeptmhmm_gff: Path
    conservation_map: Path
    signalp_output: Path
    deeptmhmm_output: Path
    comparative_output: Path


@dataclass(frozen=True)
class ExecutionResult:
    workflow: str
    written_files: dict[str, Path]
    row_counts: dict[str, int]
    classification_counts: dict[str, int]
