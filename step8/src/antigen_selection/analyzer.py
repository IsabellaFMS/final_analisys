from __future__ import annotations

from collections import defaultdict

from .config import (
    ACCESSIBLE_GROUPS,
    HIGH_STRINGENCY,
    INTERMEDIATE_STRINGENCY,
    PARTIAL_MARKERS,
    TRUNCATION_RATIO,
)
from .models import Candidate, Orthogroup, Selection, SelectionResult, TraceEntry


def build_candidates(
    accessions: set[str],
    orthogroups: dict[str, Orthogroup],
    protein_to_orthogroup: dict[str, str],
    species_by_protein: dict[str, str],
    lengths: dict[str, int],
    descriptions: dict[str, str],
    accessibility: dict[str, str],
    epitopes_high: dict[str, list[dict[str, str]]],
    epitopes_intermediate: dict[str, list[dict[str, str]]],
    mhc_class_i: dict[str, int],
    mhc_class_ii: dict[str, int],
) -> dict[str, Candidate]:
    candidates: dict[str, Candidate] = {}
    for accession in sorted(accessions):
        orthogroup = protein_to_orthogroup.get(accession, "")
        description = descriptions.get(accession, "")
        candidates[accession] = Candidate(
            accession=accession,
            description=description,
            species=species_by_protein.get(accession, ""),
            orthogroup=orthogroup,
            length=lengths.get(accession, 0),
            consensus_length=orthogroups[orthogroup].consensus_length if orthogroup in orthogroups else 0.0,
            accessibility_group=accessibility.get(accession, ""),
            epitopes_high=len(epitopes_high.get(accession, ())),
            epitopes_intermediate=len(epitopes_intermediate.get(accession, ())),
            mhc_class_i_alleles=mhc_class_i.get(accession, 0),
            mhc_class_ii_alleles=mhc_class_ii.get(accession, 0),
            partial_description=any(marker in description.lower() for marker in PARTIAL_MARKERS),
        )
    return candidates


def is_truncated(candidate: Candidate) -> bool:
    return candidate.partial_description or candidate.length_ratio < TRUNCATION_RATIO


def _rank(candidate: Candidate, epitope_count: int) -> tuple:
    # Richest epitope profile first, allelic coverage as the tie-break, then the
    # accession so the result does not depend on input order.
    return (-epitope_count, -candidate.allele_coverage, candidate.accession)


def select_antigens(
    candidates: dict[str, Candidate],
    orthogroups: dict[str, Orthogroup],
) -> SelectionResult:
    """
    Apply the five selection rules in order, recording a verdict per candidate.

    Rule 5 is not a plain exclusion: an orthogroup whose only high-stringency
    member is truncated is represented by its full-length orthologue instead,
    whose epitopes are reported at the intermediate stringency.
    """
    trace: list[TraceEntry] = []
    counts: dict[str, int] = {"proteins_considered": len(candidates)}

    def record(candidate: Candidate, verdict: str, detail: str) -> None:
        trace.append(TraceEntry(candidate.accession, candidate.orthogroup, verdict, detail))

    # Rule 1. The candidates without a high-stringency epitope stay in the
    # dictionary, because rule 5 looks for a substitute among them.
    with_evidence = [candidate for candidate in candidates.values() if candidate.epitopes_high > 0]
    counts["rule1_bcell_evidence"] = len(with_evidence)

    conserved = []
    for candidate in with_evidence:
        orthogroup = orthogroups.get(candidate.orthogroup)
        if orthogroup is None or not orthogroup.is_conserved:
            present = ", ".join(sorted(orthogroup.species)) if orthogroup else "no orthogroup"
            record(candidate, "dropped_rule2_conservation", f"species present: {present}")
        else:
            conserved.append(candidate)
    counts["rule2_conserved"] = len(conserved)

    accessible = []
    for candidate in conserved:
        if candidate.accessibility_group in ACCESSIBLE_GROUPS:
            accessible.append(candidate)
        else:
            record(candidate, "dropped_rule3_accessibility", candidate.accessibility_group or "absent from the step 5 table")
    counts["rule3_accessible"] = len(accessible)

    by_orthogroup: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in accessible:
        by_orthogroup[candidate.orthogroup].append(candidate)

    representatives: list[Candidate] = []
    for identifier, group in sorted(by_orthogroup.items()):
        ranked = sorted(group, key=lambda candidate: _rank(candidate, candidate.epitopes_high))
        representatives.append(ranked[0])
        for candidate in ranked[1:]:
            record(
                candidate,
                "dropped_rule4_representative",
                f"{ranked[0].accession} represents {identifier} with {ranked[0].epitopes_high} epitopes "
                f"and {ranked[0].allele_coverage} alleles against {candidate.epitopes_high} and {candidate.allele_coverage}",
            )
    counts["rule4_representatives"] = len(representatives)

    selected: list[Selection] = []
    substitutions = 0
    for representative in representatives:
        if not is_truncated(representative):
            record(representative, "selected", f"full length at {representative.length_ratio:.2f} of the consensus")
            selected.append(Selection(representative, HIGH_STRINGENCY, representative.epitopes_high))
            continue

        record(
            representative,
            "dropped_rule5_truncated",
            f"{representative.length} aa against a consensus of {representative.consensus_length:.0f} aa "
            f"({representative.length_ratio:.2f})",
        )
        substitute = _full_length_substitute(representative, candidates, orthogroups)
        if substitute is None:
            continue
        substitutions += 1
        record(
            substitute,
            "selected_rule5_substitute",
            f"replaces the truncated {representative.accession}; epitopes reported at {INTERMEDIATE_STRINGENCY}",
        )
        selected.append(
            Selection(substitute, INTERMEDIATE_STRINGENCY, substitute.epitopes_intermediate, replaces=representative.accession)
        )

    counts["rule5_substitutions"] = substitutions
    counts["selected"] = len(selected)
    selected.sort(key=lambda selection: selection.candidate.accession)
    return SelectionResult(
        selected=selected,
        trace=trace,
        counts=counts,
        parameters={
            "high_stringency": HIGH_STRINGENCY,
            "intermediate_stringency": INTERMEDIATE_STRINGENCY,
            "truncation_ratio": TRUNCATION_RATIO,
            "accessible_groups": "; ".join(ACCESSIBLE_GROUPS),
        },
    )


def _full_length_substitute(
    truncated: Candidate,
    candidates: dict[str, Candidate],
    orthogroups: dict[str, Orthogroup],
) -> Candidate | None:
    """
    Pick the orthologue that represents an orthogroup whose high-stringency
    member is truncated: full length, surface accessible, and carrying at least
    one epitope at the intermediate stringency.
    """
    orthogroup = orthogroups.get(truncated.orthogroup)
    if orthogroup is None:
        return None
    eligible = [
        candidate
        for member in orthogroup.members
        if (candidate := candidates.get(member)) is not None
        and candidate.accession != truncated.accession
        and not is_truncated(candidate)
        and candidate.accessibility_group in ACCESSIBLE_GROUPS
        and candidate.epitopes_intermediate > 0
    ]
    if not eligible:
        return None
    return sorted(eligible, key=lambda candidate: _rank(candidate, candidate.epitopes_intermediate))[0]
