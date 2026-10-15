from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class SpeciesConfig:
    header: str
    alias: str
    display_name: str


@dataclass
class OrthogroupRecord:
    orthogroup: str
    proteins_by_species: dict[str, list[str]]


@dataclass
class ConsolidatedRow:
    conservation_category: str
    taxonomic_group: str
    present_species_count: int
    total_orthogroups: int = 0
    total_proteins: int = 0
    strict_ortholog_proteins: int = 0
    strict_paralog_proteins: int = 0
    singleton_groups: int = 0
    singleton_proteins: int = 0
    species_orthologs: dict[str, int] = field(default_factory=dict)
    species_paralogs: dict[str, int] = field(default_factory=dict)



@dataclass(frozen=True)
class ExpandedConservationRow:
    orthogroup_id: str
    protein_id: str
    representative_protein: str
    source_virus_abbreviation: str
    compared_virus_abbreviation: str
    conservation_status: str
    viruses_with_detected_conservation: str
    viruses_without_detected_conservation: str
    orthogroup_member_count: int
    evidence_summary: str
    decision_reason: str = ""
    source_reference: str = ""


@dataclass(frozen=True)
class ProteinConservationRow:
    protein_id: str
    specie: str
    specie_name: str
    viruses_with_detected_conservation: str
    conservation: str



@dataclass
class ClassifiedProteinGroup:
    conservation_category: str
    taxonomic_group: str
    species_headers: list[str]
    present_species_count: int
    protein_ids: list[str] = field(default_factory=list)
