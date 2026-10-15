from __future__ import annotations

import logging
import unicodedata
from collections import defaultdict
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def _canonical_label(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    stripped = "".join(char for char in normalized if not unicodedata.combining(char))
    return stripped.lower()


def format_confidence(probabilities: list[str]) -> str:


    try:
        probs = [float(p) for p in probabilities if p not in {"N/A", ""}]
    except ValueError:
        return "Unknown"

    if not probs:
        return "Unknown"

    max_prob = max(probs)
    if max_prob >= 0.9:
        return "High"
    if max_prob >= 0.7:
        return "Moderate"
    return "Low"


def resolve_peptide_type(prediction: str) -> str:
    mapping = {
        "SP": "SP (Sec/SPI)",
        "LIPO": "LIPO (Sec/SPII)",
        "TAT": "TAT (Tat/SPI)",
        "TATLIPO": "TATLIPO (Tat/SPII)",
        "Pilin": "Pilin (Sec/SPIII)",
    }
    return mapping.get(prediction, prediction)


def parse_signalp(txt_path: Path, gff_path: Path) -> pd.DataFrame:
    """Parse SignalP tabular and GFF outputs into a single dataframe."""
    if not txt_path.exists() or not gff_path.exists():
        logger.warning(f"Missing SignalP inputs: {txt_path} or {gff_path}")
        return pd.DataFrame()

    predictions: dict[str, dict[str, str]] = {}
    with txt_path.open("r", encoding="utf-8") as handle:
        header: list[str] = []
        for line in handle:
            if line.startswith("#"):
                if "ID" in line and "Prediction" in line:
                    header = line.strip().replace("# ", "").split("\t")
                continue

            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            row = dict(zip(header, parts))
            pred_type = row.get("Prediction")
            if not pred_type or pred_type == "OTHER":
                continue

            prot_id = row.get("ID", "").split(" |")[0].strip()
            if not prot_id:
                continue

            sp_prob = row.get("SP(Sec/SPI)", "N/A")
            lipo_prob = row.get("LIPO(Sec/SPII)", "N/A")
            tat_prob = row.get("TAT(Tat/SPI)", "N/A")
            predictions[prot_id] = {
                "id": prot_id,
                "prediction": pred_type,
                "sp_prob": sp_prob,
                "lipo_prob": lipo_prob,
                "tat_prob": tat_prob,
                "other_prob": row.get("OTHER", "N/A"),
                "cs_position": row.get("CS Position", "N/A"),
                "n_region": "N/A",
                "h_region": "N/A",
                "c_region": "N/A",
                "confidence_level": format_confidence([sp_prob, lipo_prob, tat_prob]),
                "signal_peptide_type": resolve_peptide_type(pred_type),
            }

    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            prot_id = parts[0].split(" |")[0].strip()
            feature_type, start, end = parts[2], parts[3], parts[4]
            if prot_id not in predictions:
                continue

            coords = f"{start}-{end}"
            if feature_type == "n-region":
                predictions[prot_id]["n_region"] = coords
            elif feature_type == "h-region":
                predictions[prot_id]["h_region"] = coords
            elif feature_type == "c-region":
                predictions[prot_id]["c_region"] = coords

    return pd.DataFrame(list(predictions.values()))


def classify_deep_topology(regions: list[dict]) -> str:
    # A transmembrane helix and a signal peptide are reported first because they
    # are what decides whether the protein is accessible. Only when neither is
    # present does the label fall back to the globular forms, built from which
    # sides of the membrane the chain was placed on: E for outside, I for inside.
        
    has_tm = any(r["type"] == "TMhelix" for r in regions)
    has_sp = any(r["type"] == "signal" for r in regions)

    if has_tm and has_sp:
        return "SP+TM"
    if has_tm:
        return "TM"
    if has_sp:
        return "SP"

    has_in = any(r["type"] == "inside" for r in regions)
    has_out = any(r["type"] == "outside" for r in regions)

    parts = []
    if has_out:
        parts.append("E")
    if has_in:
        parts.append("I")
    if parts:
        parts.append("GLOB")
        return "+".join(parts)
    return "GLOB"


def parse_deeptmhmm(gff_path: Path) -> pd.DataFrame:
    """Parse DeepTMHMM GFF output into one row per protein."""
    if not gff_path.exists():
        logger.warning(f"Missing DeepTMHMM input: {gff_path}")
        return pd.DataFrame()

    region_types = ["signal", "TMhelix", "inside", "outside", "TMbeta"]
    proteins: list[dict] = []
    current: dict[str, object] = {}
    regions: list[dict[str, int | str]] = []

    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#") and "Length" in line:
                if current:
                    current["regions"] = regions
                    current["type"] = classify_deep_topology(regions)
                    proteins.append(current)
                parts = line.split()
                current = {"id": parts[1], "size": int(parts[3]), "tmrs": 0}
                regions = []
            elif line.startswith("#") and "Number of predicted TMRs" in line:
                if current:
                    current["tmrs"] = int(line.split(":")[1].strip())
            elif line.startswith("//"):
                if current:
                    current["regions"] = regions
                    current["type"] = classify_deep_topology(regions)
                    proteins.append(current)
                current, regions = {}, []
            elif not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) > 3 and parts[1] in region_types:
                    regions.append({"type": parts[1], "start": int(parts[2]), "end": int(parts[3])})

    if current:
        current["regions"] = regions
        current["type"] = classify_deep_topology(regions)
        proteins.append(current)

    rows = []
    for protein in proteins:
        row = {
            "id": protein.get("id"),
            "size": protein.get("size"),
            "tmrs": protein.get("tmrs"),
            "type": protein.get("type"),
        }
        blocks_by_type: defaultdict[str, list[str]] = defaultdict(list)
        for region in protein.get("regions", []):
            blocks_by_type[region["type"]].append(f"{region['start']}-{region['end']}")

        for region_type in region_types:
            row[f"{region_type}_present"] = "Yes" if blocks_by_type[region_type] else "No"
            for index, region in enumerate(blocks_by_type[region_type], start=1):
                row[f"{region_type}_{index}"] = region
        rows.append(row)

    return pd.DataFrame(rows)


CONSERVED_SECTION_TITLES = {"=== proteins conserved across all species ==="}
PARTIAL_SECTION_TITLES = {"=== partially conserved proteins ==="}
EXCLUSIVE_SECTION_TITLES = {"=== species-specific proteins ==="}


def build_conservation_map(txt_path: Path) -> dict[str, str]:
    """Load conservation labels from the classification report."""
    if not txt_path.exists():
        logger.warning(f"Missing conservation map: {txt_path}")
        return {}

    classifications: dict[str, str] = {}
    current_class: str | None = None
    current_partial: str | None = None
    excl_species: str | None = None

    species_aliases = {
        "cowpox_refseq_Brighton_Red_proteins": "CPXV",
        "mpox_refseq_Rivers_State_proteins": "MPXV",
        "vaccinia_WR_refseq_proteins": "VACV",
        "varv_refseq_proteins": "VARV",
    }

    with txt_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(":"):
                continue

            normalized_line = _canonical_label(line)
            if normalized_line in CONSERVED_SECTION_TITLES:
                current_class, excl_species = "Conserved", None
                continue
            if normalized_line in PARTIAL_SECTION_TITLES:
                current_class, excl_species = "Partial", None
                continue
            if normalized_line in EXCLUSIVE_SECTION_TITLES:
                current_class = "Exclusive"
                continue

            if current_class == "Partial" and ":" in line:
                species_list = line.rstrip(":").split("+")
                cleaned = [species_aliases.get(item.strip(), item.strip()) for item in species_list]
                current_partial = f'Partial ({"+".join(cleaned)})'
                continue

            if current_class == "Conserved":
                classifications[line] = "Conserved"
            elif current_class == "Partial":
                classifications[line] = current_partial or "Partial"
            elif current_class == "Exclusive":
                target_species = line.rstrip(":")
                if target_species in species_aliases:
                    excl_species = species_aliases[target_species]
                    continue
                if "_refseq_" in line:
                    excl_species = line.split("_refseq_")[0].replace("_proteins", "")
                    continue
                classifications[line] = (
                    f"Exclusive ({excl_species})" if excl_species else "Exclusive (Unknown)"
                )

    return classifications
