from __future__ import annotations

from argparse import Namespace
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CommandPaths:
    project_root: Path
    tsv_input: Path | None = None
    fasta_input: Path | None = None
    exclusion_input: Path | None = None
    exclusion_output: Path | None = None
    fasta_output: Path | None = None
    flagged_output: Path | None = None
    duplicate_report: Path | None = None
    duplicate_fasta: Path | None = None


def _resolve(project_root: Path, raw_value: str | None, default_relative: str | None) -> Path | None:
    if raw_value is None and default_relative is None:
        return None
    candidate = Path(raw_value) if raw_value is not None else project_root / default_relative
    return candidate if candidate.is_absolute() else (project_root / candidate)


def _resolve_input_with_sources_fallback(
    project_root: Path,
    raw_value: str | None,
    default_relative: str,
    glob_pattern: str,
) -> Path:
    resolved = _resolve(project_root, raw_value, default_relative)
    if raw_value is not None or resolved.exists():
        return resolved

    sources_dir = project_root / "sources"
    candidates = sorted(path for path in sources_dir.glob(glob_pattern) if path.is_file())
    if len(candidates) == 1:
        return candidates[0]
    return resolved


CLEAN_VIRAL_FASTA_NAME = "clean_viral_proteins.fasta"
DUPLICATE_REPORT_NAME = "removed_viral_proteins.csv"
DUPLICATE_FASTA_NAME = "removed_viral_proteins.fasta"


def _resolve_step3_artifact(project_root: Path, raw_value: str | None, subdir: str, name: str) -> Path:
    # The reincorporation inputs come from step 3 and are required: without them
    # the step would quietly write a proteome 45 records short.
    if raw_value is not None:
        candidate = Path(raw_value)
        return candidate if candidate.is_absolute() else (project_root / candidate)
    default = project_root.parent / "step3" / "outputs" / subdir / name
    local_copy = project_root / "sources" / name
    if not default.is_file() and local_copy.is_file():
        return local_copy
    return default


def _resolve_fasta_input_with_sources_fallback(
    project_root: Path,
    raw_value: str | None,
    default_relative: str,
) -> Path:

    resolved = _resolve(project_root, raw_value, default_relative)
    if raw_value is not None:
        return resolved

    # The proteome comes from step 3, which writes it under its own outputs
    step3_output = project_root.parent / "step3" / "outputs" / "fasta" / CLEAN_VIRAL_FASTA_NAME
    if step3_output.is_file():
        return step3_output
    if resolved.exists():
        return resolved

    sources_dir = project_root / "sources"
    candidates = sorted(
        path
        for pattern in ("*.fasta", "*.fa", "*.faa", "*.fna")
        for path in sources_dir.glob(pattern)
        if path.is_file()
    )
    if len(candidates) == 1:
        return candidates[0]

    preferred = [
        path
        for path in candidates
        if any(token in path.name.lower() for token in ("viral", "virus", "refseq"))
        and not any(token in path.name.lower() for token in ("human",))
    ]
    if len(preferred) == 1:
        return preferred[0]
    return resolved


def resolve_command_paths(command: str, project_root: Path, args: Namespace) -> CommandPaths:
    project_root = project_root.resolve()
    defaults = {
        "tsv_input": "sources/viral_hits.tsv",
        "fasta_input": f"sources/{CLEAN_VIRAL_FASTA_NAME}",
        "exclusion_input": "outputs/reports/exclusion_ids.txt",
        "exclusion_output": "outputs/reports/exclusion_ids.txt",
        "fasta_output": "outputs/fasta/safe_sequences.fasta",
        "flagged_output": "outputs/reports/flagged_reintegrated_accessions.txt",
    }
    reincorporating = command in {"run", "filter"}
    return CommandPaths(
        project_root=project_root,
        tsv_input=(
            _resolve_input_with_sources_fallback(project_root, getattr(args, "tsv_input", None), defaults["tsv_input"], "*.tsv")
            if command in {"run", "extract"}
            else None
        ),
        fasta_input=(
            _resolve_fasta_input_with_sources_fallback(project_root, getattr(args, "fasta_input", None), defaults["fasta_input"])
            if command in {"run", "filter"}
            else None
        ),
        exclusion_input=_resolve(project_root, getattr(args, "exclusion_input", None), defaults["exclusion_input"] if command == "filter" else None),
        exclusion_output=_resolve(project_root, getattr(args, "exclusion_output", None), defaults["exclusion_output"] if command in {"run", "extract"} else None),
        fasta_output=_resolve(project_root, getattr(args, "fasta_output", None), defaults["fasta_output"] if command in {"run", "filter"} else None),
        flagged_output=_resolve(project_root, getattr(args, "flagged_output", None), defaults["flagged_output"] if reincorporating else None),
        duplicate_report=(
            _resolve_step3_artifact(project_root, getattr(args, "duplicate_report", None), "reports", DUPLICATE_REPORT_NAME)
            if reincorporating
            else None
        ),
        duplicate_fasta=(
            _resolve_step3_artifact(project_root, getattr(args, "duplicate_fasta", None), "fasta", DUPLICATE_FASTA_NAME)
            if reincorporating
            else None
        ),
    )
