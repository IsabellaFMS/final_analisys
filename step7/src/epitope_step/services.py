from __future__ import annotations

from .analyzer import (
    add_protein_ids,
    analyze_mhcii_promiscuity,
    build_mhci_final_output,
    build_mhcii_final_output,
    concatenate_tables,
    filter_hits_by_rank,
)
from .config import StepPaths
from .loaders import load_fasta_id_map, load_summary_classification_map, read_table
from .models import WorkflowResult
from .writer import write_csv


def _log(message: str) -> None:
    print(message, flush=True)


def run_workflows(
    *,
    paths: StepPaths,
    workflows: tuple[str, ...],
    mhci_cutoff: float,
    mhcii_cutoffs: tuple[float, ...],
) -> list[WorkflowResult]:
    results: list[WorkflowResult] = []
    for workflow in workflows:
        if workflow == "mhci":
            _log("[mhci] Starting workflow.")
            results.append(run_mhci_pipeline(paths, mhci_cutoff))
            _log("[mhci] Workflow finished.")
        elif workflow == "mhcii":
            _log("[mhcii] Starting workflow.")
            results.append(run_mhcii_pipeline(paths, mhcii_cutoffs))
            _log("[mhcii] Workflow finished.")
        else:
            raise ValueError(f"Unsupported workflow: {workflow}")
    return results


def run_mhci_pipeline(paths: StepPaths, cutoff: float) -> WorkflowResult:
    input_files = sorted(paths.mhci_prediction_dir.glob("prediction_*.csv"))
    if not input_files:
        raise FileNotFoundError(f"No MHCI prediction files found in {paths.mhci_prediction_dir}")
    _log(f"[mhci] Reading {len(input_files)} prediction files from sources.")
    frames = [read_table(path, sep=r"\s+") for path in input_files]
    _log(f"[mhci] Filtering hits with rank <= {cutoff:g}.")
    filtered = filter_hits_by_rank(concatenate_tables(frames), cutoff)
    _log(f"[mhci] Mapping protein IDs from FASTA and building final table ({len(filtered)} filtered hits).")
    ids = add_protein_ids(filtered, load_fasta_id_map(paths.fasta_file))
    classification_map = load_summary_classification_map(paths.classification_file)
    final_output = build_mhci_final_output(ids, classification_map)
    artifact = write_csv(final_output, paths.mhci_final_output)
    _log(f"[mhci] Wrote final output: {artifact.name} ({len(final_output)} rows).")

    return WorkflowResult(
        workflow="mhci",
        artifacts=[artifact],
        metrics={"input_files": len(input_files), "final_rows": len(final_output)},
        final_frame=final_output,
    )


def run_mhcii_pipeline(paths: StepPaths, cutoffs: tuple[float, ...]) -> WorkflowResult:


    if not paths.mhcii_prediction_files:
        raise FileNotFoundError("No MHCII prediction files configured.")
    _log(f"[mhcii] Reading {len(paths.mhcii_prediction_files)} prediction files from sources.")
    base_frame = concatenate_tables([read_table(path) for path in paths.mhcii_prediction_files])
    id_map = load_fasta_id_map(paths.fasta_file)
    enriched_frames = []

    for cutoff in cutoffs:
        _log(f"[mhcii] Applying cutoff {cutoff:g} and collecting enriched hits.")
        filtered = filter_hits_by_rank(base_frame, cutoff)
        enriched_frames.append(add_protein_ids(filtered, id_map))

    _log("[mhcii] Consolidating cutoff results and annotating allele counts per protein/epitope.")
    merged = concatenate_tables(enriched_frames, deduplicate=True)
    promiscuous = analyze_mhcii_promiscuity(merged)
    classification_map = load_summary_classification_map(paths.classification_file)
    final_output = build_mhcii_final_output(promiscuous, classification_map)
    artifact = write_csv(final_output, paths.mhcii_final_output)
    _log(f"[mhcii] Wrote final output: {artifact.name} ({len(final_output)} rows).")

    return WorkflowResult(
        workflow="mhcii",
        artifacts=[artifact],
        metrics={"cutoffs": len(cutoffs), "final_rows": len(final_output)},
        final_frame=final_output,
    )
