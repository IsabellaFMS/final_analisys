import re
from pathlib import Path

import pandas as pd

from .config import SUMMARY_GO_COLUMNS, SUMMARY_REACTOME_COLUMNS
from .writer import write_text


def build_basic_report(summary_df: pd.DataFrame, stats: dict[str, int | float]) -> str:
    lines = [
        "DETAILED PROTEOME ANALYSIS",
        f"Total Reference Proteins: {stats['total_proteins']}",
        "",
        "1. GENERAL QUANTITATIVE DATA (Criterion: GO term presence)",
        f"GO annotation success rate: {stats['annotated_proteins']/stats['total_proteins']:.2%} ({stats['annotated_proteins']} of {stats['total_proteins']} proteins).",
        f"Total proteins WITHOUT GO annotation in any tool: {stats['unannotated_proteins']}",
        "-" * 60,
        "Overall InterProScan vs. eggNOG comparison (based on GO terms):",
        f"- With GO from eggNOG: {stats['eggnog_go']}",
        f"- With GO from InterProScan: {stats['interpro_go']}",
        f"- With GO in BOTH: {stats['both_go']}",
        f"- GAIN | GO exclusive to eggNOG: {stats['eggnog_only_go']}",
        f"- GAIN | GO exclusive to InterProScan: {stats['interpro_only_go']}",
        "-" * 60,
        "",
        "3-7. DETAILED SUMMARY TABLE BY CONSERVATION GROUP",
        "",
        summary_df[["Group", "Total Proteins", "With GO (%)", "Top 5 Domains (InterPro)", "Top 5 GOs (Combined)"]].to_string(index=False),
    ]
    return "\n".join(lines) + "\n"


def build_enhanced_report(summary_df: pd.DataFrame, total_proteins: int) -> str:
    return "\n".join([
        "ENHANCED PROTEOME ANALYSIS",
        f"Total Reference Proteins: {total_proteins}",
        "",
        "DETAILED SUMMARY TABLE BY CONSERVATION GROUP",
        "",
        summary_df.to_string(index=False),
    ]) + "\n"


def annotate_summary_ids(df: pd.DataFrame, go_map: dict[str, str], reactome_map: dict[str, str]) -> pd.DataFrame:
    annotated = df.copy()

    def replace_ids(text: str, resolver: dict[str, str]) -> str:
        if not isinstance(text, str) or text == "None found":
            return text
        return re.sub(r"([^,()]+)\s*\((\d+)\)", lambda m: f"{resolver.get(m.group(1).strip(), m.group(1).strip())} ({m.group(2)})", text)

    for column in SUMMARY_GO_COLUMNS:
        if column in annotated.columns:
            annotated[column] = annotated[column].apply(lambda value: replace_ids(value, go_map))
    for column in SUMMARY_REACTOME_COLUMNS:
        if column in annotated.columns:
            annotated[column] = annotated[column].apply(lambda value: replace_ids(value, reactome_map))
    return annotated


def replace_ranked_items(text: str, resolver: dict[str, str]) -> str:
    if not isinstance(text, str) or text == "None found":
        return text
    items = []
    for item in text.split(","):
        stripped = item.strip()
        match = re.fullmatch(r"(.+?)\s*\((\d+)\)", stripped)
        if not match:
            items.append(stripped)
            continue
        label = match.group(1).strip()
        count = match.group(2)
        items.append(f"{resolver.get(label, label)} ({count})")
    return ", ".join(items)


def resolve_summary_labels(df: pd.DataFrame, go_map: dict[str, str] | None = None) -> pd.DataFrame:
    resolved = df.copy()
    go_map = go_map or {}

    for column in SUMMARY_GO_COLUMNS:
        if column in resolved.columns:
            resolved[column] = resolved[column].apply(lambda value: replace_ranked_items(value, go_map))

    return resolved


def write_enrichment_report(enrichment_files: list[Path], output_path: Path) -> Path:
    lines = ["CONSOLIDATED GO TERM ENRICHMENT REPORT", ""]
    for file_path in enrichment_files:
        df = pd.read_csv(file_path, sep="\t")
        df["p_fdr_bh"] = pd.to_numeric(df["p_fdr_bh"])
        group_name = file_path.stem.replace("enrichment_", "")
        lines.append(f"Group: {group_name}")
        lines.append(f"Total significant terms: {len(df)}")
        lines.append(f"Enriched: {(df['enrichment'] == 'e').sum()}")
        lines.append(f"Purified: {(df['enrichment'] == 'p').sum()}")
        lines.append("")
    return write_text("\n".join(lines), output_path)
