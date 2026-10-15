from __future__ import annotations

from orthogroup_analysis.models import ClassifiedProteinGroup, ConsolidatedRow, ExpandedConservationRow, OrthogroupRecord, ProteinConservationRow, SpeciesConfig



CATEGORY_ORDER = {
    "Conserved": 0,
    "Partially_Conserved": 1,
    "Species-Specific": 2,
}


def analyze_orthogroups(
    records: list[OrthogroupRecord],
    species_config: list[SpeciesConfig],
) -> list[ConsolidatedRow]:
    alias_by_header = {species.header: species.alias for species in species_config}
    total_species = len(species_config)
    consolidated: dict[str, ConsolidatedRow] = {}


    for record in records:
        present_headers = [
            header for header, proteins in record.proteins_by_species.items() if proteins
        ]


        if not present_headers:
            continue

        present_aliases = [alias_by_header[header] for header in present_headers]
        present_species_count = len(present_aliases)
        taxonomic_group = " + ".join(present_aliases)

        conservation_category = classify_conservation(present_species_count, total_species)

        row = consolidated.get(taxonomic_group)
        if row is None:
            row = ConsolidatedRow(
                conservation_category=conservation_category,
                taxonomic_group=taxonomic_group,
                present_species_count=present_species_count,
              
                species_orthologs={species.alias: 0 for species in species_config},
                species_paralogs={species.alias: 0 for species in species_config},
            )
            consolidated[taxonomic_group] = row

        total_proteins = 0
        strict_orthologs = 0
        strict_paralogs = 0

        for header in present_headers:
            alias = alias_by_header[header]
            protein_count = len(record.proteins_by_species[header])

            total_proteins += protein_count

            # One strict ortholog per species present, regardless of how many copies 
            # that species contributes. The extra copies are counted as paralogs
            strict_orthologs += 1
            strict_paralogs += max(0, protein_count - 1)

            row.species_orthologs[alias] += 1
            row.species_paralogs[alias] += max(0, protein_count - 1)

        row.total_orthogroups += 1
        row.total_proteins += total_proteins
        row.strict_ortholog_proteins += strict_orthologs
        row.strict_paralog_proteins += strict_paralogs

        if conservation_category == "Species-Specific":
            row.singleton_groups += 1
            row.singleton_proteins += total_proteins

    return sorted(consolidated.values(), key=sort_key)


def build_expanded_conservation_rows(
    records: list[OrthogroupRecord],
    species_config: list[SpeciesConfig],
) -> list[ExpandedConservationRow]:


    alias_by_header = {species.header: species.alias for species in species_config}
    all_aliases = [species.alias for species in species_config]
    rows: list[ExpandedConservationRow] = []

    for record in records:
        present_headers = [
            header for header, proteins in record.proteins_by_species.items() if proteins
        ]
        if not present_headers:
            continue

        present_aliases = [alias_by_header[header] for header in present_headers]
        absent_aliases = [alias for alias in all_aliases if alias not in present_aliases]
        conservation_category = classify_conservation(
            len(present_aliases),
            len(species_config),
        )
        representative_protein = select_representative_protein(record, species_config)
        member_count = count_members(record, present_headers)
        viruses_with_conservation = ";".join(present_aliases)
        viruses_without_conservation = ";".join(absent_aliases)

        for species in species_config:
            source_proteins = record.proteins_by_species.get(species.header, [])
            if not source_proteins:
                continue

            for protein_id in source_proteins:
                for compared_alias in all_aliases:
                    if compared_alias == species.alias:
                        continue

                    pair_is_conserved = compared_alias in present_aliases
                    rows.append(
                        ExpandedConservationRow(
                            orthogroup_id=record.orthogroup,
                            protein_id=protein_id,
                            representative_protein=representative_protein,
                            source_virus_abbreviation=species.alias,
                            compared_virus_abbreviation=compared_alias,
                            conservation_status=conservation_category,
                            viruses_with_detected_conservation=viruses_with_conservation,
                            viruses_without_detected_conservation=viruses_without_conservation,
                            orthogroup_member_count=member_count,
                            evidence_summary=build_evidence_summary(
                                record.orthogroup,
                                protein_id,
                                species.alias,
                                compared_alias,
                                present_aliases,
                                pair_is_conserved,
                            ),
                            decision_reason=build_decision_reason(
                                species.alias,
                                compared_alias,
                                pair_is_conserved,
                            ),
                            source_reference="step1_orthogroup_presence_analysis",
                        )
                    )

    return rows


def build_protein_conservation_rows(
    records: list[OrthogroupRecord],
    species_config: list[SpeciesConfig],
) -> list[ProteinConservationRow]:
    species_by_header = {species.header: species for species in species_config}
    rows: list[ProteinConservationRow] = []

    for record in records:
        present_headers = [
            header for header, proteins in record.proteins_by_species.items() if proteins
        ]
        if not present_headers:
            continue

        present_aliases = [species_by_header[header].alias for header in present_headers]
        conservation_category = classify_conservation(
            len(present_headers),
            len(species_config),
        )
        viruses_with_detected_conservation = ";".join(present_aliases)

        for header in present_headers:
            species = species_by_header[header]
            for protein_id in record.proteins_by_species[header]:
                rows.append(
                    ProteinConservationRow(
                        protein_id=protein_id,
                        specie=species.alias,
                        specie_name=species.display_name,
                        viruses_with_detected_conservation=viruses_with_detected_conservation,
                        conservation=conservation_category,
                    )
                )

    return rows


def build_classified_protein_groups(
    records: list[OrthogroupRecord],
    species_config: list[SpeciesConfig],
) -> list[ClassifiedProteinGroup]:
    alias_by_header = {species.header: species.alias for species in species_config}
    groups: dict[str, ClassifiedProteinGroup] = {}

    for record in records:
        present_headers = [
            header for header, proteins in record.proteins_by_species.items() if proteins
        ]
        if not present_headers:
            continue

        ordered_headers = [
            species.header
            for species in species_config
            if species.header in present_headers
        ]
        present_aliases = [alias_by_header[header] for header in ordered_headers]
        taxonomic_group = " + ".join(present_aliases)

        group = groups.get(taxonomic_group)
        if group is None:
            group = ClassifiedProteinGroup(
                conservation_category=classify_conservation(
                    len(ordered_headers),
                    len(species_config),
                ),
                taxonomic_group=taxonomic_group,
                species_headers=ordered_headers,
                present_species_count=len(ordered_headers),
            )
            groups[taxonomic_group] = group

        for header in ordered_headers:
            group.protein_ids.extend(record.proteins_by_species[header])

    return sorted(groups.values(), key=classified_group_sort_key)


def classified_group_sort_key(group: ClassifiedProteinGroup) -> tuple[int, int, str]:
    return (
        CATEGORY_ORDER[group.conservation_category],
        -group.present_species_count,
        group.taxonomic_group,
    )


def classify_conservation(present_species_count: int, total_species: int) -> str:
    if present_species_count == total_species:
        return "Conserved"

    if present_species_count == 1:
        return "Species-Specific"

    return "Partially_Conserved"


def sort_key(row: ConsolidatedRow) -> tuple[int, int, str]:
    return (
        CATEGORY_ORDER[row.conservation_category],
        -row.present_species_count,
        row.taxonomic_group,
    )


def select_representative_protein(
    record: OrthogroupRecord,
    species_config: list[SpeciesConfig],
) -> str:

   
    for species in species_config:
        proteins = record.proteins_by_species.get(species.header, [])
        if proteins:
            return proteins[0]
    return ""


def count_members(record: OrthogroupRecord, present_headers: list[str]) -> int:
    return sum(len(record.proteins_by_species[header]) for header in present_headers)


def build_evidence_summary(
    orthogroup_id: str,
    protein_id: str,
    source_alias: str,
    compared_alias: str,
    present_aliases: list[str],
    pair_is_conserved: bool,
) -> str:
    if pair_is_conserved:
        return (
            f"Protein {protein_id} belongs to orthogroup {orthogroup_id}, which "
            f"contains detected members for {', '.join(present_aliases)}; therefore "
            f"the {source_alias} to {compared_alias} pair is conserved by orthogroup presence."
        )

    return (
        f"Protein {protein_id} belongs to orthogroup {orthogroup_id}, which contains "
        f"detected members for {', '.join(present_aliases)} but none for {compared_alias}; "
        f"therefore the {source_alias} to {compared_alias} pair is not conserved."
    )


def build_decision_reason(
    source_alias: str,
    compared_alias: str,
    pair_is_conserved: bool,
) -> str:
    if pair_is_conserved:
        return (
            f"The orthogroup contains at least one {compared_alias} member alongside "
            f"the {source_alias} source protein."
        )

    return (
        f"The orthogroup does not contain any {compared_alias} member to match "
        f"the {source_alias} source protein."
    )
