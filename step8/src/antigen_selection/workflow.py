from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from .analyzer import build_candidates, select_antigens
from .config import HIGH_STRINGENCY, INTERMEDIATE_STRINGENCY, sibling_step
from .models import SelectionResult, StepPaths
from .loaders import (
    load_accessibility,
    load_allele_coverage,
    load_epitopes,
    load_orthogroups,
    load_sequences,
    require,
)
from .writer import write_antigens, write_selected_epitopes, write_trace


def epitope_table(paths: StepPaths, cutoff: float) -> Path:
    return paths.bcell / f"{cutoff}" / "tables" / "epitopes.csv"


def extract_epitopes(paths: StepPaths, cutoff: float, refresh: bool = False) -> Path:
    """
    Run the step 6 extraction at one cutoff, into this step's output tree.

    Step 6 is deposited at 0.1512 and its table must not be overwritten, so the
    runs this step needs are directed elsewhere through its own options.
    """
    table = epitope_table(paths, cutoff)
    if table.is_file() and not refresh:
        return table

    step6 = sibling_step(paths, "step6")
    require(step6 / "run_analysis.py", "Step 6 runner", "step 6")
    subprocess.run(
        [
            sys.executable, "run_analysis.py",
            "--cutoff", str(cutoff),
            "--output-dir", str(paths.bcell / f"{cutoff}"),
        ],
        cwd=step6,
        check=True,
    )
    return require(table, f"Step 6 epitope table at cutoff {cutoff}", "step 6")


def run_selection(paths: StepPaths, refresh: bool = False) -> tuple[SelectionResult, dict[float, dict[str, list[dict[str, str]]]]]:
    epitopes = {
        cutoff: load_epitopes(extract_epitopes(paths, cutoff, refresh))
        for cutoff in (HIGH_STRINGENCY, INTERMEDIATE_STRINGENCY)
    }

    accessibility = load_accessibility(
        require(
            sibling_step(paths, "step5", "outputs", "tables", "comparative_results_signalp6_deep_final.csv"),
            "Step 5 comparative table",
            "step 5",
        )
    )
    lengths, descriptions = load_sequences(
        require(
            sibling_step(paths, "step2", "sources", "all_refseq_proteins.fasta"),
            "Reference proteome",
            "step 2",
        )
    )
    orthogroups, protein_to_orthogroup, species_by_protein = load_orthogroups(
        require(sibling_step(paths, "step1", "sources", "orthogroups.tsv"), "Orthogroup table", "step 1"),
        lengths,
    )
    mhc_class_i = load_allele_coverage(
        require(sibling_step(paths, "step7", "outputs", "mhci_final.csv"), "MHC-I table", "step 7")
    )
    mhc_class_ii = load_allele_coverage(
        require(sibling_step(paths, "step7", "outputs", "mhcii_final.csv"), "MHC-II table", "step 7")
    )

    candidates = build_candidates(
        accessions=set(accessibility) | {accession for table in epitopes.values() for accession in table},
        orthogroups=orthogroups,
        protein_to_orthogroup=protein_to_orthogroup,
        species_by_protein=species_by_protein,
        lengths=lengths,
        descriptions=descriptions,
        accessibility=accessibility,
        epitopes_high=epitopes[HIGH_STRINGENCY],
        epitopes_intermediate=epitopes[INTERMEDIATE_STRINGENCY],
        mhc_class_i=mhc_class_i,
        mhc_class_ii=mhc_class_ii,
    )

    result = select_antigens(candidates, orthogroups)
    write_antigens(result, paths.tables / "final_antigens.tsv")
    write_trace(result, paths.tables / "selection_trace.tsv")
    write_selected_epitopes(result.selected, epitopes, paths.tables / "selected_epitopes.tsv")
    return result, epitopes
