from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


CLASSIFICATION_FILE_NAME = "orthogroups_consolidated_analysis_classified_protein_ids.txt"


def resolve_classification_file(root: Path) -> Path:

    root = root.resolve()
    step1_report = root.parent / "step1" / "outputs" / CLASSIFICATION_FILE_NAME
    if step1_report.is_file():
        return step1_report
    return root / "sources" / CLASSIFICATION_FILE_NAME


def format_cutoff(cutoff: float) -> str:
    return f"{cutoff:g}"


@dataclass(frozen=True)
class StepPaths:
    root: Path
    sources: Path
    outputs: Path
    docs: Path
    tests: Path
    fasta_file: Path
    classification_file: Path
    mhci_prediction_dir: Path
    mhcii_prediction_files: tuple[Path, ...]

    @classmethod
    def from_root(cls, root: Path) -> "StepPaths":
        root = root.resolve()
        sources = root / "sources"
        return cls(
            root=root,
            sources=sources,
            outputs=root / "outputs",
            docs=root / "docs",
            tests=root / "tests",
            fasta_file=sources / "safe_sequences.fasta",
            classification_file=resolve_classification_file(root),
            mhci_prediction_dir=sources / "mhci" / "full_predictions",
            mhcii_prediction_files=(
                sources / "mhcii" / "DRB1" / "DRB1_predict_result.csv",
                sources / "mhcii" / "DQA1" / "DQA1_predict_result.csv",
            ),
        )

    @property
    def mhci_final_output(self) -> Path:
        return self.outputs / "mhci_final.csv"

    @property
    def mhcii_final_output(self) -> Path:
        return self.outputs / "mhcii_final.csv"
