from __future__ import annotations

import re
import unicodedata


ALL_SPECIES = "CPXV, MPXV, VACV, VARV"
SPECIES_ABBREVIATIONS = {
    "cowpox_refseq_Brighton_Red_proteins": "CPXV",
    "mpox_refseq_Rivers_State_proteins": "MPXV",
    "varv_refseq_proteins": "VARV",
    "vaccinia_WR_refseq_proteins": "VACV",
}


def canonical_text(text: str) -> str:
    normalized = unicodedata.normalize("NFKD", text)
    stripped = "".join(char for char in normalized if not unicodedata.combining(char))
    return stripped.lower().strip()


def translate_conservation_label(label: str) -> str:
    normalized = canonical_text(label)
    if normalized == "unclassified":
        return "Unclassified"

    if normalized in {
        "proteins conserved across all species",
        "conserved"
    }:
        return f"Conserved ({ALL_SPECIES})"

    if (
        normalized.startswith("partially conserved proteins:") or
        normalized.startswith("partially conserved:")
    ):
        content = label.split(":", 1)[1]
        species = _extract_species_list(content)
        return f"Partially conserved ({species})"

    if (
        normalized.startswith("species-specific proteins:") or
        normalized.startswith("species-specific:")
    ):
        content = label.split(":", 1)[1]
        species = _extract_species_list(content)
        return f"Exclusive ({species})"

    if normalized.startswith("exclusive "):
        content = label.split(" ", 1)[1]
        species = _extract_species_list(content)
        return f"Exclusive ({species})"

    if "+" in label:
        return f"Partially conserved ({_extract_species_list(label)})"
    return label


def normalize_species_name(raw_species: str | None) -> str:
    if not raw_species:
        return "Unknown"
    species = raw_species.strip()
    if species.lower().endswith(" virus"):
        return f"Orthopoxvirus {species[:-6].strip().lower()}"
    return species


def extract_species_from_header(header: str) -> str | None:
    match = re.search(r"\[([^\]]+)\]\s*$", header)
    if not match:
        return None
    return normalize_species_name(match.group(1))


def _extract_species_list(raw_text: str) -> str:
    cleaned = raw_text.replace(":", "")
    for long_name, short_name in SPECIES_ABBREVIATIONS.items():
        cleaned = cleaned.replace(long_name, short_name)
    tokens = [token.strip() for token in cleaned.split("+") if token.strip()]
    return ", ".join(tokens)
