from __future__ import annotations

import csv
from pathlib import Path

from orthogroup_analysis.models import ClassifiedProteinGroup, ConsolidatedRow, ExpandedConservationRow, ProteinConservationRow, SpeciesConfig



CATEGORY_SECTION_TITLES = {
    "Conserved": "Proteins Conserved Across All Species",
    "Partially_Conserved": "Partially Conserved Proteins",
    "Species-Specific": "Species-Specific Proteins",
}


def build_fieldnames(species_config: list[SpeciesConfig]) -> list[str]:
    fieldnames = [
        "conservationCategory",
        "taxonomicGroup",
        "presentSpeciesCount",
        "totalOrthogroups",
        "totalProteins",
        "strictOrthologProteins",
        "strictParalogProteins",
        "singletonGroups",
        "singletonProteins",
    ]

    for species in species_config:
        fieldnames.append(f"{species.alias}_Orthologs")
        fieldnames.append(f"{species.alias}_Paralogs")

    return fieldnames


def build_row_payloads(
    rows: list[ConsolidatedRow],
    species_config: list[SpeciesConfig],
) -> list[dict[str, str | int]]:
    payloads: list[dict[str, str | int]] = []
    for row in rows:
        payload = {
            "conservationCategory": row.conservation_category,
            "taxonomicGroup": row.taxonomic_group,
            "presentSpeciesCount": row.present_species_count,
            "totalOrthogroups": row.total_orthogroups,
            "totalProteins": row.total_proteins,
            "strictOrthologProteins": row.strict_ortholog_proteins,
            "strictParalogProteins": row.strict_paralog_proteins,
            "singletonGroups": row.singleton_groups,
            "singletonProteins": row.singleton_proteins,
        }

        for species in species_config:
            payload[f"{species.alias}_Orthologs"] = row.species_orthologs[species.alias]
            payload[f"{species.alias}_Paralogs"] = row.species_paralogs[species.alias]

        payloads.append(payload)

    return payloads


def write_consolidated_csv(
    output_path: Path,
    rows: list[ConsolidatedRow],
    species_config: list[SpeciesConfig],
) -> None:
    fieldnames = build_fieldnames(species_config)
    payloads = build_row_payloads(rows, species_config)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for payload in payloads:
            writer.writerow(payload)


def build_expanded_fieldnames() -> list[str]:
    return [
        "orthogroup_id",
        "protein_id",
        "representative_protein",
        "source_virus_abbreviation",
        "compared_virus_abbreviation",
        "conservation_status",
        "viruses_with_detected_conservation",
        "viruses_without_detected_conservation",
        "orthogroup_member_count",
        "evidence_summary",
        "decision_reason",
        "source_reference",
    ]


def build_expanded_row_payloads(
    rows: list[ExpandedConservationRow],
) -> list[dict[str, str | int]]:
    return [
        {
            "orthogroup_id": row.orthogroup_id,
            "protein_id": row.protein_id,
            "representative_protein": row.representative_protein,
            "source_virus_abbreviation": row.source_virus_abbreviation,
            "compared_virus_abbreviation": row.compared_virus_abbreviation,
            "conservation_status": row.conservation_status,
            "viruses_with_detected_conservation": row.viruses_with_detected_conservation,
            "viruses_without_detected_conservation": row.viruses_without_detected_conservation,
            "orthogroup_member_count": row.orthogroup_member_count,
            "evidence_summary": row.evidence_summary,
            "decision_reason": row.decision_reason,
            "source_reference": row.source_reference,
        }
        for row in rows
    ]


def write_expanded_conservation_csv(
    output_path: Path,
    rows: list[ExpandedConservationRow],
) -> None:
    fieldnames = build_expanded_fieldnames()
    payloads = build_expanded_row_payloads(rows)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for payload in payloads:
            writer.writerow(payload)


def build_protein_conservation_fieldnames() -> list[str]:
    return [
        "protein_id",
        "specie",
        "specie_name",
        "viruses_with_detected_conservation",
        "conservation",
    ]


def build_protein_conservation_payloads(
    rows: list[ProteinConservationRow],
) -> list[dict[str, str]]:
    return [
        {
            "protein_id": row.protein_id,
            "specie": row.specie,
            "specie_name": row.specie_name,
            "viruses_with_detected_conservation": row.viruses_with_detected_conservation,
            "conservation": row.conservation,
        }
        for row in rows
    ]


def write_protein_conservation_csv(
    output_path: Path,
    rows: list[ProteinConservationRow],
) -> None:
    fieldnames = build_protein_conservation_fieldnames()
    payloads = build_protein_conservation_payloads(rows)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for payload in payloads:
            writer.writerow(payload)


def build_classified_protein_ids_lines(groups: list[ClassifiedProteinGroup]) -> list[str]:
    lines: list[str] = []
    current_category: str | None = None

    for group in groups:
        if group.conservation_category != current_category:
            current_category = group.conservation_category

            if lines:
                lines.append("")

            lines.append(f"=== {CATEGORY_SECTION_TITLES[current_category]} ===")

        if current_category != "Conserved":
            lines.append("")
            lines.append(" + ".join(group.species_headers) + ":")

        lines.extend(group.protein_ids)

    return lines


def write_classified_protein_ids_txt(
    output_path: Path,
    groups: list[ClassifiedProteinGroup],
) -> None:
    lines = build_classified_protein_ids_lines(groups)

    output_path.parent.mkdir(parents=True, exist_ok=True)


    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for line in lines:
            handle.write(f"{line}\n")
