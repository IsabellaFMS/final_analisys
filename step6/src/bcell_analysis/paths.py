from __future__ import annotations

from pathlib import Path

from .models import StepPaths


def _resolve_path(base_dir: Path, candidate: Path) -> Path:
    return candidate.resolve() if candidate.is_absolute() else (base_dir / candidate).resolve()


CLASSIFICATION_FILE_NAME = "orthogroups_consolidated_analysis_classified_protein_ids.txt"


def resolve_classification_file(project_root: Path) -> Path:
    # Step 1 is the only producer; the copy under sources/ is the fallback.
    step1_report = project_root.parent / "step1" / "outputs" / CLASSIFICATION_FILE_NAME
    if step1_report.is_file():
        return step1_report
    return project_root / "sources" / CLASSIFICATION_FILE_NAME


def build_step_paths(
    project_root: Path,
    input_file: Path,
    output_dir: Path,
    classification_file: Path | None = None,
) -> StepPaths:
    root = project_root.resolve()
    resolved_input = _resolve_path(root, input_file)
    resolved_output = _resolve_path(root, output_dir)
    conservation_map = (
        _resolve_path(root, classification_file)
        if classification_file is not None
        else resolve_classification_file(root)
    )
    return StepPaths(
        project_root=root,
        input_file=resolved_input,
        output_root=resolved_output,
        tables_dir=resolved_output / "tables",
        fasta_dir=resolved_output / "fasta",
        conservation_map=conservation_map,
    )
