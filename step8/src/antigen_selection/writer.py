from __future__ import annotations

import csv
from pathlib import Path

from .models import Selection, SelectionResult

ANTIGEN_COLUMNS = [
    "Protein_ID", "Description", "Species", "Orthogroup", "Length",
    "Consensus_Length", "Length_Ratio", "Reported_Stringency", "Epitope_Count",
    "MHC_I_Alleles", "MHC_II_Alleles", "Accessibility_Group", "Replaces",
]

TRACE_COLUMNS = ["Protein_ID", "Orthogroup", "Verdict", "Detail"]


def _prepare_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _write_rows(path: Path, columns: list[str], rows: list[dict[str, object]]) -> Path:
    _prepare_parent(path)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_antigens(result: SelectionResult, path: Path) -> Path:
    rows = []
    for selection in result.selected:
        candidate = selection.candidate
        rows.append({
            "Protein_ID": candidate.accession,
            "Description": candidate.description,
            "Species": candidate.species,
            "Orthogroup": candidate.orthogroup,
            "Length": candidate.length,
            "Consensus_Length": f"{candidate.consensus_length:.0f}",
            "Length_Ratio": f"{candidate.length_ratio:.2f}",
            "Reported_Stringency": selection.stringency,
            "Epitope_Count": selection.epitope_count,
            "MHC_I_Alleles": candidate.mhc_class_i_alleles,
            "MHC_II_Alleles": candidate.mhc_class_ii_alleles,
            "Accessibility_Group": candidate.accessibility_group,
            "Replaces": selection.replaces or "-",
        })
    return _write_rows(path, ANTIGEN_COLUMNS, rows)


def write_trace(result: SelectionResult, path: Path) -> Path:
    rows = [
        {"Protein_ID": entry.accession, "Orthogroup": entry.orthogroup, "Verdict": entry.verdict, "Detail": entry.detail}
        for entry in sorted(result.trace, key=lambda entry: (entry.orthogroup, entry.accession))
    ]
    return _write_rows(path, TRACE_COLUMNS, rows)


def write_selected_epitopes(
    selections: list[Selection],
    epitopes_by_stringency: dict[float, dict[str, list[dict[str, str]]]],
    path: Path,
) -> Path:
    rows: list[dict[str, object]] = []
    columns: list[str] = []
    for selection in selections:
        for epitope in epitopes_by_stringency[selection.stringency].get(selection.candidate.accession, ()):
            if not columns:
                columns = list(epitope.keys()) + ["Cutoff"]
            rows.append({**epitope, "Cutoff": selection.stringency})
    return _write_rows(path, columns or ["Cutoff"], rows)
