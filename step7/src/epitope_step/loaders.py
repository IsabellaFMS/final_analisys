from __future__ import annotations

from pathlib import Path

import pandas as pd

from .translation import canonical_text, extract_species_from_header, translate_conservation_label


ABBREVIATIONS = {
    "cowpox_refseq_Brighton_Red_proteins": "CPXV",
    "mpox_refseq_Rivers_State_proteins": "MPXV",
    "varv_refseq_proteins": "VARV",
    "vaccinia_WR_refseq_proteins": "VACV",
}


def read_table(path: Path, *, sep: str | None = None) -> pd.DataFrame:
    if sep is None:
        return pd.read_csv(path)
    return pd.read_csv(path, sep=sep)


def load_fasta_id_map(path: Path) -> dict[int, str]:
    mapping: dict[int, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        index = 0
        for line in handle:
            if line.startswith(">"):
                index += 1
                mapping[index] = line[1:].strip().split()[0]
    return mapping


def load_fasta_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            parts = header.split(None, 1)
            protein_id = parts[0]
            description = parts[1] if len(parts) > 1 else ""
            metadata[protein_id] = {
                "description": description,
                "species": extract_species_from_header(header) or "Unknown",
            }
    return metadata


def load_summary_classification_map(path: Path) -> dict[str, str]:
    classification_map: dict[str, str] = {}
    current_classification = "Unclassified"
    main_category = ""

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("==="):
                main_category = line.strip("= ").strip()
                current_classification = translate_conservation_label(main_category)
                continue
            if "+" in line or (":" in line and not line.startswith(("NP_", "YP_"))):
                sub_category = line.replace(":", "").strip()
                current_classification = translate_conservation_label(f"{main_category}: {sub_category}")
                continue
            protein_id = line.split()[0]
            classification_map[protein_id] = current_classification
    return classification_map


def load_hit_classification_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    main_class: str | None = None
    sub_class: str | None = None

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("==="):
                main_class = line.strip("= ").strip()
                sub_class = None
                continue
            normalized_main = canonical_text(main_class or "")
            if normalized_main in {
                "partially conserved proteins",
                "species-specific proteins",
            } and ("+" in line or line.endswith(":")):
                label = line.strip(":").strip()
                for long_name, short_name in ABBREVIATIONS.items():
                    label = label.replace(long_name, short_name)
                sub_class = " + ".join(part.strip() for part in label.split("+"))
                continue
            final_class = translate_conservation_label(main_class or "Unclassified")
            if sub_class:
                if normalized_main == "species-specific proteins":
                    final_class = translate_conservation_label(f"Exclusive {sub_class}")
                else:
                    final_class = translate_conservation_label(sub_class)
            mapping[line.split()[0]] = final_class
    return mapping
