from __future__ import annotations

from .models import SelectionResult

RULE_LABELS = {
    "proteins_considered": "Proteins considered (step 4 proteome)",
    "rule1_bcell_evidence": "Rule 1, at least one high-stringency B-cell epitope",
    "rule2_conserved": "Rule 2, orthogroup conserved across the four species",
    "rule3_accessible": "Rule 3, surface accessible by consensus",
    "rule4_representatives": "Rule 4, one representative per orthogroup",
    "rule5_substitutions": "Rule 5, truncated representatives replaced",
    "selected": "Selected antigens",
}


def format_summary(result: SelectionResult) -> str:
    lines = ["Final antigen selection", ""]
    lines.append("Parameters")
    for key, value in result.parameters.items():
        lines.append(f"  {key}: {value}")

    lines.extend(["", "Rules"])
    for key, label in RULE_LABELS.items():
        if key in result.counts:
            lines.append(f"  {label}: {result.counts[key]}")

    lines.extend(["", "Selected"])
    for selection in result.selected:
        candidate = selection.candidate
        lines.append(
            f"  {candidate.accession}  {candidate.description or '-'} "
            f"[{candidate.species}, {candidate.orthogroup}]"
        )
        lines.append(
            f"    {selection.epitope_count} epitopes at {selection.stringency}, "
            f"{candidate.mhc_class_i_alleles} MHC-I and {candidate.mhc_class_ii_alleles} MHC-II alleles"
        )
        if selection.replaces:
            lines.append(f"    represents its orthogroup in place of the truncated {selection.replaces}")

    truncated = [entry for entry in result.trace if entry.verdict == "dropped_rule5_truncated"]
    if truncated:
        lines.extend(["", "Rule 5 fired for"])
        for entry in truncated:
            lines.append(f"  {entry.orthogroup}: {entry.accession}, {entry.detail}")

    return "\n".join(lines) + "\n"
