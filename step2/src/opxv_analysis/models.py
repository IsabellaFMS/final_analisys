from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PipelinePaths:

    root: Path
    sources: Path
    outputs: Path
    tables: Path
    reports: Path
    texts: Path
    visualizations: Path
    goatools: Path
    enrichment: Path


@dataclass(frozen=True)
class ComparisonStats:

    eggnog_unique: int
    interpro_unique: int
    union_total: int
    both: int
    interpro_only: int
    eggnog_only: int
