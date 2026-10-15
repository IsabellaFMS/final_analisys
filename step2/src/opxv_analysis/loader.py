import re
from pathlib import Path

import pandas as pd
from goatools.obo_parser import GODag

from .config import CONSERVATION_SPECIES, EGGNOG_COLUMNS
from .models import PipelinePaths
from .writer import ensure_dirs


def setup_paths(paths: PipelinePaths) -> PipelinePaths:
    ensure_dirs(paths.outputs, paths.tables, paths.reports, paths.texts, paths.visualizations, paths.goatools, paths.enrichment)
    return paths


def load_eggnog(paths: PipelinePaths) -> pd.DataFrame:
    return pd.read_csv(paths.sources / "eggnog_annotations.tsv", sep="\t", comment="#", header=None, names=EGGNOG_COLUMNS, dtype=str)


def load_interpro(paths: PipelinePaths) -> pd.DataFrame:
    return pd.read_csv(paths.sources / "interpro_orthopoxvirus.tsv", sep="\t", dtype=str)


def parse_conservation_list_name(stem: str) -> tuple[str, list[str]] | None:
    # Map the txt file to its conservation group and species
    if stem == "list_core":
        return "core", ["Undetermined"]
    for prefix in ("list_exclusive_", "list_partial_"):
        if stem.startswith(prefix):
            species = sorted(part for part in stem[len(prefix):].split("_") if part in CONSERVATION_SPECIES)
            if species:
                return "+".join(species), species
    return None


def load_conservation_groups(paths: PipelinePaths) -> pd.DataFrame:
    
    records = []
    for path in sorted(paths.sources.glob("list_*.txt")):
        parsed = parse_conservation_list_name(path.stem)
        if parsed is None:
            continue
        group, species = parsed
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            protein_id = line.strip()
            if not protein_id:
                continue
            for name in species:
                records.append({"Protein_ID": protein_id, "Conservation_Group": group, "Species": name})
    if not records:
        raise FileNotFoundError(f"No list_*.txt conservation file was found in {paths.sources}")
    return pd.DataFrame(records).drop_duplicates()


def load_go_dag(paths: PipelinePaths) -> GODag:
    return GODag(str(paths.sources / "go-basic.obo"), optional_attrs={"relationship"})


def load_go_name_map(paths: PipelinePaths) -> dict[str, str]:
    # Build a compact GO ID translation to term name dictionary from the OBO file
    go_map: dict[str, str] = {}
    current_id = None
    for line in (paths.sources / "go-basic.obo").read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("id:"):
            current_id = line.split(" ", 1)[1].strip()
        elif line.startswith("name:") and current_id:
            go_map[current_id] = line.split(" ", 1)[1].strip()
            current_id = None
    return go_map


def load_reactome_map(paths: PipelinePaths) -> dict[str, str]:
    candidates = [paths.sources / "ReactomePathways.txt", paths.root / "ReactomePathways.txt"]
    for candidate in candidates:
        if candidate.exists():
            df = pd.read_csv(candidate, sep="\t", header=None, usecols=[0, 1], names=["id", "name"])
            return pd.Series(df.name.values, index=df.id).to_dict()
    return {}


def load_fasta_count(paths: PipelinePaths) -> int:
    with open(paths.sources / "all_refseq_proteins.fasta", "r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def find_enrichment_files(search_root: Path) -> list[Path]:

    return sorted(search_root.glob("enrichment_*.tsv"))


def sanitize_group_name(group_name: str) -> str:
    return group_name.replace("+", "_plus_")


def extract_go_ids(value: str) -> list[str]:
    return re.findall(r"GO:\d{7}", str(value))
