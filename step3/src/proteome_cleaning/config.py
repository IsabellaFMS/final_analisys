from __future__ import annotations

from argparse import Namespace
from dataclasses import dataclass
from pathlib import Path


DEFAULT_KMER_SIZE = 12

DEFAULTS = {
    "human_input": "sources/human.fasta",
    "viral_input": "sources/viral_proteins.fasta",
    "human_output": "outputs/fasta/clean_human_proteome.fasta",
    "viral_output": "outputs/fasta/clean_viral_proteins.fasta",
    "human_report": "outputs/reports/removed_human_proteins.csv",
    "viral_report": "outputs/reports/removed_viral_proteins.csv",
    "viral_removed_output": "outputs/fasta/removed_viral_proteins.fasta",
    "summary_output": "outputs/reports/cleaning_summary.txt",
}

HUMAN_COMMANDS = {"run", "clean-human"}
VIRAL_COMMANDS = {"run", "clean-viral"}


@dataclass(frozen=True)
class CommandPaths:
    project_root: Path
    kmer_size: int = DEFAULT_KMER_SIZE
    human_input: Path | None = None
    viral_input: Path | None = None
    human_output: Path | None = None
    viral_output: Path | None = None
    human_report: Path | None = None
    viral_report: Path | None = None
    viral_removed_output: Path | None = None
    summary_output: Path | None = None


def _resolve(project_root: Path, raw_value: str | None, default_relative: str | None) -> Path | None:
    if raw_value is None and default_relative is None:
        return None
    candidate = Path(raw_value) if raw_value is not None else project_root / default_relative
    return candidate if candidate.is_absolute() else (project_root / candidate)


def _default_for(command: str, key: str) -> str | None:
    wanted = HUMAN_COMMANDS if key.startswith("human") else VIRAL_COMMANDS
    return DEFAULTS[key] if command in wanted else None


def resolve_command_paths(command: str, project_root: Path, args: Namespace) -> CommandPaths:
    project_root = project_root.resolve()
    kmer_size = getattr(args, "kmer_size", None) or DEFAULT_KMER_SIZE
    if kmer_size < 1:
        raise ValueError(f"K-mer size must be at least 1, got {kmer_size}.")
    keys = (
        "human_input", "viral_input", "human_output", "viral_output",
        "human_report", "viral_report", "viral_removed_output",
    )
    resolved = {
        key: _resolve(project_root, getattr(args, key, None), _default_for(command, key))
        for key in keys
    }
    return CommandPaths(
        project_root=project_root,
        kmer_size=kmer_size,
        summary_output=_resolve(
            project_root,
            getattr(args, "summary_output", None),
            DEFAULTS["summary_output"] if command == "run" else None,
        ),
        **resolved,
    )
