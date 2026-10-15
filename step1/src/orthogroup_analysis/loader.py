from __future__ import annotations

import csv
from pathlib import Path

from orthogroup_analysis.config import ORTHOGROUP_COLUMN, build_species_config
from orthogroup_analysis.models import OrthogroupRecord, SpeciesConfig


def load_records(input_paths: list[Path]) -> tuple[list[SpeciesConfig], list[OrthogroupRecord]]:
    species_config: list[SpeciesConfig] | None = None

    records_by_orthogroup: dict[str, OrthogroupRecord] = {}

    for path in input_paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            headers = reader.fieldnames or []

            validate_headers(path, headers)

            current_species_config = build_species_config(headers)
            if species_config is None:
                species_config = current_species_config
            else:
                validate_species_compatibility(path, species_config, current_species_config)

            for row in reader:
                proteins_by_species = {
                    species.header: parse_protein_cell(row.get(species.header, ""))
                    for species in current_species_config
                }

                orthogroup = (row.get(ORTHOGROUP_COLUMN) or "").strip()
                merge_record(
                    records_by_orthogroup,
                    OrthogroupRecord(
                        orthogroup=orthogroup,
                        proteins_by_species=proteins_by_species,
                    ),
                )

    if species_config is None:
        raise ValueError("No input files were loaded.")

    records = sorted(records_by_orthogroup.values(), key=lambda item: item.orthogroup)
    return species_config, records


def validate_headers(path: Path, headers: list[str]) -> None:
    if ORTHOGROUP_COLUMN not in headers:
        raise ValueError(f"File {path} is missing required column '{ORTHOGROUP_COLUMN}'.")

    if len(headers) < 2:
        raise ValueError(f"File {path} does not contain any species columns.")


def validate_species_compatibility(
    path: Path,
    expected: list[SpeciesConfig],
    current: list[SpeciesConfig],
) -> None:
    expected_headers = [item.header for item in expected]
    current_headers = [item.header for item in current]
    if expected_headers != current_headers:
        raise ValueError(
            f"File {path} has species columns different from the first input file: "
            f"{current_headers} != {expected_headers}"
        )


def parse_protein_cell(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def merge_record(
    records_by_orthogroup: dict[str, OrthogroupRecord],
    incoming: OrthogroupRecord,
) -> None:
    current = records_by_orthogroup.get(incoming.orthogroup)
    if current is None:
        records_by_orthogroup[incoming.orthogroup] = incoming
        return

    for species_header, proteins in incoming.proteins_by_species.items():
        merged = current.proteins_by_species[species_header]
        for protein in proteins:
            if protein not in merged:
                merged.append(protein)
