from __future__ import annotations

from pathlib import Path

from .models import StepPaths

DEFAULT_SIGNALP_TXT = Path("sources/signalp6/prediction_results.txt")
DEFAULT_SIGNALP_GFF = Path("sources/signalp6/region_output.gff3")
DEFAULT_DEEP_GFF = Path("sources/deeptmhmm/TMRs.gff3")

CLASSIFICATION_FILE_NAME = "orthogroups_consolidated_analysis_classified_protein_ids.txt"

# Columns kept in the comparative table, in order, mapped to their output name.
COMPARATIVE_COLUMNS = {
    "id": "id",
    "prediction": "Prediction SignalP",
    "type": "Prediction Deep",
    "strategic_classification": "strategic_classification",
    "classification": "classification",
}

DEFAULT_SIGNALP_OUT = Path("outputs/tables/signalp6_results_final.csv")
DEFAULT_DEEP_OUT = Path("outputs/tables/deeptmhmm_summary_final.csv")
DEFAULT_FINAL_MERGE = Path("outputs/tables/comparative_results_signalp6_deep_final.csv")


def resolve_classification_file(project_root: Path) -> Path:
    project_root = project_root.resolve()
    step1_report = project_root.parent / "step1" / "outputs" / CLASSIFICATION_FILE_NAME
    if step1_report.is_file():
        return step1_report
    return project_root / "sources" / CLASSIFICATION_FILE_NAME


def _resolve_from_root(project_root: Path, candidate: Path | None, default: Path) -> Path:
    target = candidate if candidate is not None else default
    return target.resolve() if target.is_absolute() else (project_root / target).resolve()


def resolve_step_paths(
    project_root: Path,
    workflow: str = "full",
    signalp_txt: Path | None = None,
    signalp_gff: Path | None = None,
    deeptmhmm_gff: Path | None = None,
    conservation_map: Path | None = None,
    signalp_output: Path | None = None,
    deeptmhmm_output: Path | None = None,
    comparative_output: Path | None = None,
) -> StepPaths:
    project_root = project_root.resolve()
    return StepPaths(
        project_root=project_root,
        workflow=workflow,
        sources_dir=(project_root / "sources").resolve(),
        outputs_dir=(project_root / "outputs").resolve(),
        docs_dir=(project_root / "docs").resolve(),
        signalp_txt=_resolve_from_root(project_root, signalp_txt, DEFAULT_SIGNALP_TXT),
        signalp_gff=_resolve_from_root(project_root, signalp_gff, DEFAULT_SIGNALP_GFF),
        deeptmhmm_gff=_resolve_from_root(project_root, deeptmhmm_gff, DEFAULT_DEEP_GFF),
        conservation_map=_resolve_from_root(
            project_root,
            conservation_map,
            resolve_classification_file(project_root),
        ),
        signalp_output=_resolve_from_root(project_root, signalp_output, DEFAULT_SIGNALP_OUT),
        deeptmhmm_output=_resolve_from_root(project_root, deeptmhmm_output, DEFAULT_DEEP_OUT),
        comparative_output=_resolve_from_root(project_root, comparative_output, DEFAULT_FINAL_MERGE),
    )
