from __future__ import annotations

import logging

import pandas as pd

from .analyzer import annotate_dataframe, merge_analytics
from .config import COMPARATIVE_COLUMNS
from .loaders import build_conservation_map, parse_deeptmhmm, parse_signalp
from .models import ExecutionResult, StepPaths
from .writer import export_csv


def _configure_logging() -> None:
    logging.basicConfig(
        level=logging.WARNING,
        format="%(levelname)s %(name)s %(message)s",
        handlers=[logging.StreamHandler()],
        force=True,
    )


def _classification_counts(df: pd.DataFrame) -> dict[str, int]:
    if "classification" not in df.columns:
        return {}

    counts = df["classification"].fillna("Unclassified").value_counts(dropna=False)
    return {str(label): int(count) for label, count in counts.items()}


def _signalp_table(paths: StepPaths) -> pd.DataFrame:
    signal_df = parse_signalp(paths.signalp_txt, paths.signalp_gff)
    cons_map = build_conservation_map(paths.conservation_map)
    return annotate_dataframe(signal_df, cons_map)


def _deeptmhmm_table(paths: StepPaths) -> pd.DataFrame:
    deep_df = parse_deeptmhmm(paths.deeptmhmm_gff)
    cons_map = build_conservation_map(paths.conservation_map)
    return annotate_dataframe(deep_df, cons_map)


def _comparative_table(paths: StepPaths) -> pd.DataFrame:
    signal_df = parse_signalp(paths.signalp_txt, paths.signalp_gff)
    deep_df = parse_deeptmhmm(paths.deeptmhmm_gff)
    cons_map = build_conservation_map(paths.conservation_map)
    merged = annotate_dataframe(merge_analytics(signal_df, deep_df), cons_map)
    return _project_comparative(merged)


def _project_comparative(df: pd.DataFrame) -> pd.DataFrame:
    # confidence_level and signal_present stay out of this table, but they are
    # still computed: signal_present is what strategic_classification is built
    # from, and both remain in the SignalP and DeepTMHMM tables.
    if df.empty:
        return pd.DataFrame(columns=list(COMPARATIVE_COLUMNS.values()))
    return df.reindex(columns=list(COMPARATIVE_COLUMNS)).rename(columns=COMPARATIVE_COLUMNS)


def run_signalp_workflow(paths: StepPaths) -> ExecutionResult:
    _configure_logging()
    signal_df = _signalp_table(paths)
    export_csv(signal_df, paths.signalp_output)
    return ExecutionResult(
        workflow="signalp",
        written_files={"signalp": paths.signalp_output},
        row_counts={"signalp": len(signal_df)},
        classification_counts=_classification_counts(signal_df),
    )


def run_deeptmhmm_workflow(paths: StepPaths) -> ExecutionResult:
    _configure_logging()
    deep_df = _deeptmhmm_table(paths)
    export_csv(deep_df, paths.deeptmhmm_output)
    return ExecutionResult(
        workflow="deeptmhmm",
        written_files={"deeptmhmm": paths.deeptmhmm_output},
        row_counts={"deeptmhmm": len(deep_df)},
        classification_counts=_classification_counts(deep_df),
    )


def run_comparative_workflow(paths: StepPaths) -> ExecutionResult:
    _configure_logging()
    comparative_df = _comparative_table(paths)
    export_csv(comparative_df, paths.comparative_output)
    return ExecutionResult(
        workflow="comparative",
        written_files={"comparative": paths.comparative_output},
        row_counts={"comparative": len(comparative_df)},
        classification_counts=_classification_counts(comparative_df),
    )


def run_full_workflow(paths: StepPaths) -> ExecutionResult:
    _configure_logging()
    signal_df = _signalp_table(paths)
    deep_df = _deeptmhmm_table(paths)
    comparative_df = _comparative_table(paths)

    export_csv(signal_df, paths.signalp_output)
    export_csv(deep_df, paths.deeptmhmm_output)
    export_csv(comparative_df, paths.comparative_output)

    return ExecutionResult(
        workflow="full",
        written_files={
            "signalp": paths.signalp_output,
            "deeptmhmm": paths.deeptmhmm_output,
            "comparative": paths.comparative_output,
        },
        row_counts={
            "signalp": len(signal_df),
            "deeptmhmm": len(deep_df),
            "comparative": len(comparative_df),
        },
        classification_counts=_classification_counts(comparative_df),
    )
