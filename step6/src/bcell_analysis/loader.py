from __future__ import annotations

import csv
import unicodedata
from pathlib import Path
from typing import Dict, Iterable, List

from .models import ProteinRecord

ACCESSION_FIELD = "Accession"
RESIDUE_FIELD = "Residue"
SCORE_FIELD = "BepiPred-3.0 linear epitope score"


def _clean_row(row: dict[str, str | None]) -> dict[str, str]:
    cleaned = {
        (key or "").strip(): (value or "").strip()
        for key, value in row.items()
    }
    cleaned.pop("", None)
    return cleaned


def _validate_headers(fieldnames: Iterable[str] | None) -> None:
    if not fieldnames:
        raise ValueError("Input file is empty or missing a header row.")

    available = {name.strip() for name in fieldnames if name}
    required = {ACCESSION_FIELD, RESIDUE_FIELD, SCORE_FIELD}
    missing = sorted(required - available)
    if missing:
        raise ValueError(f"Input file is missing required columns: {', '.join(missing)}")


def load_bepipred_data(input_file: Path, separator: str) -> List[ProteinRecord]:
    """Group the residue-level BepiPred-3.0 export into one record per protein."""
    proteins: Dict[str, ProteinRecord] = {}

    with input_file.open("r", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter=separator)
        _validate_headers(reader.fieldnames)

        for row in reader:
            record = _clean_row(row)
            accession = record.get(ACCESSION_FIELD, "")
            if not accession:
                continue

            protein = proteins.setdefault(accession, ProteinRecord(accession=accession))
            protein.residues += record.get(RESIDUE_FIELD, "")

            try:
                protein.scores.append(float(record.get(SCORE_FIELD, "0")))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid numeric score for accession {accession!r}."
                ) from exc

    return list(proteins.values())


CONSERVED_SECTION_TITLE = "=== proteins conserved across all species ==="
PARTIAL_SECTION_TITLE = "=== partially conserved proteins ==="
EXCLUSIVE_SECTION_TITLE = "=== species-specific proteins ==="

SPECIES_ALIASES = {
    "cowpox_refseq_Brighton_Red_proteins": "CPXV",
    "mpox_refseq_Rivers_State_proteins": "MPXV",
    "vaccinia_WR_refseq_proteins": "VACV",
    "varv_refseq_proteins": "VARV",
}

UNCLASSIFIED_LABEL = "Unclassified"


def _canonical_label(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    return "".join(char for char in normalized if not unicodedata.combining(char)).lower()


def build_conservation_map(txt_path: Path) -> Dict[str, str]:
    # The step 1 report, read the same way step 5 and step 7 read it: section
    # headers set the class, and the accessions below them inherit it.
    if not txt_path.exists():
        raise FileNotFoundError(
            f"Conservation report not found: {txt_path}. It is written by step 1, "
            "so run step 1 first, or pass --classification-file."
        )

    classifications: Dict[str, str] = {}
    current_class: str | None = None
    current_partial: str | None = None
    exclusive_species: str | None = None

    with txt_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith(":"):
                continue

            normalized = _canonical_label(line)
            if normalized == CONSERVED_SECTION_TITLE:
                current_class, exclusive_species = "Conserved", None
                continue
            if normalized == PARTIAL_SECTION_TITLE:
                current_class, exclusive_species = "Partial", None
                continue
            if normalized == EXCLUSIVE_SECTION_TITLE:
                current_class = "Exclusive"
                continue

            if current_class == "Partial" and ":" in line:
                species = [SPECIES_ALIASES.get(item.strip(), item.strip()) for item in line.rstrip(":").split("+")]
                current_partial = f'Partial ({"+".join(species)})'
                continue

            if current_class == "Conserved":
                classifications[line] = "Conserved"
            elif current_class == "Partial":
                classifications[line] = current_partial or "Partial"
            elif current_class == "Exclusive":
                target = line.rstrip(":")
                if target in SPECIES_ALIASES:
                    exclusive_species = SPECIES_ALIASES[target]
                    continue
                if "_refseq_" in line:
                    exclusive_species = line.split("_refseq_")[0].replace("_proteins", "")
                    continue
                classifications[line] = f"Exclusive ({exclusive_species})" if exclusive_species else "Exclusive (Unknown)"

    if not classifications:
        raise ValueError(f"No accession could be read from the conservation report: {txt_path}")
    return classifications


def protein_accession(value: str) -> str:
    # The BepiPred export keeps the full FASTA header; the report keys are bare.
    return value.split(" |")[0].split("|")[0].strip()
