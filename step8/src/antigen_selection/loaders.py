from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

from .config import OUTGROUP_COLUMNS, SPECIES_BY_COLUMN
from .models import Orthogroup


def _accession(value: str) -> str:
    # The step 6 and step 7 tables carry the full FASTA header in their
    # identifier column, so the accession is its first field.
    return str(value).replace("|", " ").split()[0] if str(value).strip() else ""


def require(path: Path, label: str, step: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(
            f"{label} not found: {path}. It is written by {step}, so run {step} "
            "first, or point at the file with the matching command-line option."
        )
    return path


def load_epitopes(path: Path) -> dict[str, list[dict[str, str]]]:
    """Group the step 6 epitope table by protein."""
    epitopes: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter=";"):
            epitopes[_accession(row["Accession"])].append(row)
    return dict(epitopes)


def load_accessibility(path: Path) -> dict[str, str]:
    """Read the step 5 comparative table as protein to strategic classification."""
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["id"].strip(): row["strategic_classification"].strip() for row in csv.DictReader(handle, delimiter=";")}


def load_allele_coverage(path: Path) -> dict[str, int]:
    """Count the distinct alleles binding any peptide of each protein."""
    alleles: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = _accession(row["Protein_ID"])
            for allele in str(row.get("Alleles", "")).split(","):
                if allele.strip():
                    alleles[accession].add(allele.strip())
    return {accession: len(values) for accession, values in alleles.items()}


def load_sequences(path: Path) -> tuple[dict[str, int], dict[str, str]]:
    """Read the proteome once, returning sequence lengths and descriptions."""
    lengths: dict[str, int] = {}
    descriptions: dict[str, str] = {}
    accession = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            accession = header.split()[0]
            lengths[accession] = 0
            descriptions[accession] = header[len(accession):].strip().lstrip("|").strip()
        elif accession is not None:
            lengths[accession] += len(line.strip())
    return lengths, descriptions


def load_orthogroups(path: Path, lengths: dict[str, int]) -> tuple[dict[str, Orthogroup], dict[str, str], dict[str, str]]:
    """
    Read the orthogroup table, dropping the outgroup columns.

    The consensus length is the median length of the members, so a single
    truncated or extended member does not move it.
    """
    orthogroups: dict[str, Orthogroup] = {}
    protein_to_orthogroup: dict[str, str] = {}
    species_by_protein: dict[str, str] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        keep = [index for index, name in enumerate(header) if index > 0 and name.strip() not in OUTGROUP_COLUMNS]
        for row in reader:
            if not row:
                continue
            identifier = row[0].strip()
            members: list[str] = []
            present: set[str] = set()
            for index in keep:
                if index >= len(row):
                    continue
                accessions = [value.strip() for value in row[index].split(",") if value.strip()]
                if not accessions:
                    continue
                species = SPECIES_BY_COLUMN.get(header[index].strip(), header[index].strip())
                present.add(species)
                for accession in accessions:
                    members.append(accession)
                    protein_to_orthogroup[accession] = identifier
                    species_by_protein[accession] = species
            if not members:
                continue
            member_lengths = [lengths[member] for member in members if member in lengths]
            orthogroups[identifier] = Orthogroup(
                identifier=identifier,
                members=tuple(members),
                species=frozenset(present),
                consensus_length=statistics.median(member_lengths) if member_lengths else 0.0,
            )
    return orthogroups, protein_to_orthogroup, species_by_protein
