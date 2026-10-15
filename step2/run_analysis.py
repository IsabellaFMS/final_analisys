from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import subprocess
import sys

import pandas as pd

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from opxv_analysis.analyzer import (
    build_goatools_inputs,
    build_group_summary,
    compare_annotation_coverage,
    consolidate_annotations,
    merge_annotations,
    refine_annotations,
)
from opxv_analysis.config import default_paths
from opxv_analysis.enrichment import run_dual_level_enrichment
from opxv_analysis.loader import (
    find_enrichment_files,
    load_conservation_groups,
    load_eggnog,
    load_fasta_count,
    load_go_dag,
    load_go_name_map,
    load_interpro,
    load_reactome_map,
    setup_paths,
)
from opxv_analysis.reporting import annotate_summary_ids, build_basic_report, build_enhanced_report, resolve_summary_labels, write_enrichment_report
from opxv_analysis.writer import write_table, write_text


def resolve_rscript() -> str:
    from_path = shutil.which("Rscript")
    if from_path:
        return from_path

    configured = os.environ.get("OPXV_RSCRIPT")
    if configured and Path(configured).exists():
        return configured

    raise SystemExit(
        "Rscript was not found. Install R and add Rscript to PATH, or set "
        "OPXV_RSCRIPT to its full path, before using the heatmap or full commands."
    )


def run_heatmap() -> None:
    script_path = Path(__file__).resolve().with_name("orthopoxvirus_go_heatmap.R")
    if not script_path.exists():
        raise SystemExit(f"Heatmap script not found: {script_path}")

    paths = default_paths()
    (paths.visualizations / "orthopoxvirus_go_heatmap.html").unlink(missing_ok=True)

    print("Running the R heatmap step...", flush=True)
    command = [resolve_rscript(), str(script_path)]
    env = os.environ.copy()
    env["OPXV_PROJECT_ROOT"] = str(default_paths().root)
    env["LANGUAGE"] = "en"  # keeps R warnings/errors in English on localized systems
    subprocess.run(command, check=True, env=env)


def build_summary_dataframe(paths) -> tuple[pd.DataFrame, dict[str, int | float]]:
    eggnog = load_eggnog(paths)
    interpro = load_interpro(paths)
    conservation_df = load_conservation_groups(paths)
    return build_group_summary(eggnog, interpro, conservation_df, load_fasta_count(paths))


def run_pipeline() -> None:
    paths = setup_paths(default_paths())
    eggnog = load_eggnog(paths)
    interpro = load_interpro(paths)
    go_dag = load_go_dag(paths)
    reactome_map = load_reactome_map(paths)

    conservation_df = load_conservation_groups(paths)
    (paths.tables / "detailed_conservation_groups.tsv").unlink(missing_ok=True)

    # The merged and refined tables are intermediate steps towards the
    # consolidated table, so only the latter is written to disk.
    (paths.tables / "step1_eggnog_interpro_merge.tsv").unlink(missing_ok=True)
    (paths.tables / "step2_refined_annotations.tsv").unlink(missing_ok=True)

    merged = merge_annotations(eggnog, interpro)
    refined = refine_annotations(merged, go_dag, reactome_map)
    consolidated = consolidate_annotations(refined, go_dag)
    write_table(consolidated, paths.tables / "step3_final_consolidated_annotations.tsv")

    summary_df, stats = build_group_summary(eggnog, interpro, conservation_df, load_fasta_count(paths))
    labeled_summary_df = resolve_summary_labels(summary_df, load_go_name_map(paths))
    write_table(labeled_summary_df, paths.tables / "enhanced_summary_table.tsv")
    (paths.tables / "translated_enhanced_summary_table.tsv").unlink(missing_ok=True)
    write_text(build_basic_report(summary_df, stats), paths.reports / "final_report.txt")
    write_text(build_enhanced_report(summary_df, stats["total_proteins"]), paths.reports / "enhanced_final_report.txt")
    build_goatools_inputs(eggnog, interpro, conservation_df, paths)
    run_dual_level_enrichment(paths)
    print("Python analysis completed.", flush=True)


def run_full_workflow() -> None:
    print("Running the complete workflow: Python analysis followed by the R heatmap step...", flush=True)
    run_pipeline()
    run_heatmap()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the OPXV analysis workflow.")
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("pipeline")
    sub.add_parser("heatmap")
    sub.add_parser("full")
    sub.add_parser("compare")
    sub.add_parser("annotate-summary")
    sub.add_parser("enrichment-report")
    sub.add_parser("enrichment")

    args = parser.parse_args()
    command = args.command or "full"

    paths = setup_paths(default_paths())
    if command == "pipeline":
        run_pipeline()
    elif command == "heatmap":
        run_heatmap()
    elif command == "full":
        run_full_workflow()
    elif command == "compare":
        stats = compare_annotation_coverage(load_eggnog(paths), load_interpro(paths))
        print(stats)
    elif command == "annotate-summary":
        summary, _ = build_summary_dataframe(paths)
        annotated = annotate_summary_ids(summary, load_go_name_map(paths), load_reactome_map(paths))
        write_table(annotated, paths.tables / "named_summary_table.tsv")
    elif command == "enrichment-report":
        files = find_enrichment_files(paths.root)
        if not files:
            raise SystemExit("No enrichment_*.tsv file was found.")
        write_enrichment_report(files, paths.reports / "final_enrichment_report.txt")
    elif command == "enrichment":
        run_dual_level_enrichment(paths)
    else:
        parser.error(f"Unknown command: {command}")


if __name__ == "__main__":
    main()
