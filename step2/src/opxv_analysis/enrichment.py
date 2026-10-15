from __future__ import annotations

import copy
import csv
import sys
import types
from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable

from .config import OUTGROUP_COLUMNS
from .models import PipelinePaths


def _install_multipletests_shim() -> None:
    # goatools imports multipletests from statsmodels.sandbox.stats.multicomp.
    # Recent statsmodels keeps that module but no longer exports the name, so the
    # attribute has to be checked rather than the import, and installed before
    # goatools is imported.
    try:
        import statsmodels.stats.multitest as multitest
    except ImportError:
        return
    try:
        import statsmodels.sandbox.stats.multicomp as sandbox
    except ImportError:
        sandbox = types.ModuleType("statsmodels.sandbox.stats.multicomp")
        sys.modules["statsmodels.sandbox.stats.multicomp"] = sandbox
    if not hasattr(sandbox, "multipletests"):
        sandbox.multipletests = multitest.multipletests


_install_multipletests_shim()

from goatools.go_enrichment import GOEnrichmentStudy  # noqa: E402
from goatools.obo_parser import GODag  # noqa: E402

ALPHA = 0.05
CORE_GROUP = "core"
AGGREGATION_RULES = ("union", "intersection")

RESULT_COLUMNS = [
    "GO", "name", "namespace", "enrichment", "study_count", "study_n",
    "pop_count", "pop_n", "fold_enrichment", "p_uncorrected", "p_fdr_bh",
]


def load_population(path: Path) -> dict[str, set[str]]:
    """Read goatools_population.txt, one '<accession>\\t<GO:x;GO:y>' per line."""
    associations: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        accession, _, terms = line.partition("\t")
        accession = accession.strip()
        go_terms = {term.strip() for term in terms.split(";") if term.strip()}
        if accession and go_terms:
            associations[accession] = go_terms
    return associations


def load_study_files(directory: Path) -> dict[str, list[str]]:
    """Map each conservation group to the accessions of its study file."""
    groups: dict[str, list[str]] = {}
    for path in sorted(directory.glob("goatools_study_*.txt")):
        group = path.stem[len("goatools_study_"):]
        groups[group] = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return groups


def load_protein_to_orthogroup(path: Path) -> dict[str, str]:
    """Read the protein-to-orthogroup map, dropping the outgroup columns."""
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        keep = [index for index, name in enumerate(header) if index > 0 and name.strip() not in OUTGROUP_COLUMNS]
        for row in reader:
            if not row:
                continue
            orthogroup = row[0].strip()
            for index in keep:
                if index >= len(row):
                    continue
                for accession in row[index].split(","):
                    accession = accession.strip()
                    if accession:
                        mapping[accession] = orthogroup
    return mapping


def aggregate(
    associations: dict[str, set[str]],
    protein_to_orthogroup: dict[str, str],
    rule: str,
) -> tuple[dict[str, set[str]], dict[str, list[str]]]:
    """
    Collapse a protein-level association map into an orthogroup-level one.

    Only annotated proteins contribute, so `intersection` is the intersection
    across the annotated members rather than across every member.
    """
    if rule not in AGGREGATION_RULES:
        raise ValueError(f"Unknown aggregation rule: {rule}")

    members: dict[str, list[str]] = defaultdict(list)
    for accession in associations:
        orthogroup = protein_to_orthogroup.get(accession)
        if orthogroup is not None:
            members[orthogroup].append(accession)

    aggregated: dict[str, set[str]] = {}
    for orthogroup, accessions in members.items():
        term_sets = [associations[accession] for accession in accessions]
        if rule == "union":
            terms = set().union(*term_sets)
        else:
            terms = set(term_sets[0]).intersection(*term_sets[1:])
        if terms:
            aggregated[orthogroup] = terms
    return aggregated, members


def _fold_enrichment(study: tuple[int, int], population: tuple[int, int]) -> float:
    try:
        return (study[0] / study[1]) / (population[0] / population[1])
    except ZeroDivisionError:
        return float("nan")


def run_enrichment(
    study_ids: Iterable[str],
    associations: dict[str, set[str]],
    dag: GODag,
    output_path: Path,
    log: Callable[[str], None],
    alpha: float = ALPHA,
) -> dict[str, object]:
    """Run one enrichment and write every tested term, not only the significant ones."""
    population = set(associations)
    study_ids = list(study_ids)
    # The study set has to be restricted to the population before the test:
    # passing ids that carry no annotation inflates the study size and distorts
    # the p-values.
    study = sorted(set(study_ids) & population)
    dropped = len(set(study_ids) - population)

    log(f"    population={len(population)}  study_in={len(study_ids)}  study_used={len(study)}  dropped_not_in_population={dropped}")

    summary: dict[str, object] = {
        "population": len(population),
        "study_in": len(study_ids),
        "study_used": len(study),
        "dropped": dropped,
        "terms_tested": 0,
        "significant": 0,
        "top_term": "",
        "top_name": "",
        "top_fe": "",
        "top_p": "",
        "top_fdr": "",
    }
    if not study:
        log("    study is empty after restriction to the population; skipped")
        return summary

    # propagate_counts=True mutates the association dict in place, so GOATOOLS
    # gets a copy and the caller keeps a pristine one.
    working = copy.deepcopy({key: set(value) for key, value in associations.items()})
    study_object = GOEnrichmentStudy(
        list(working.keys()),
        working,
        dag,
        propagate_counts=True,
        alpha=alpha,
        methods=["fdr_bh"],
    )

    rows = []
    for result in study_object.run_study(study):
        in_study = getattr(result, "ratio_in_study", (0, 0))
        in_population = getattr(result, "ratio_in_pop", (0, 0))
        rows.append({
            "GO": result.GO,
            "name": result.name,
            "namespace": getattr(result, "NS", ""),
            "enrichment": getattr(result, "enrichment", ""),
            "study_count": in_study[0],
            "study_n": in_study[1],
            "pop_count": in_population[0],
            "pop_n": in_population[1],
            "fold_enrichment": f"{_fold_enrichment(in_study, in_population):.6g}",
            "p_uncorrected": f"{result.p_uncorrected:.6e}",
            "p_fdr_bh": f"{result.p_fdr_bh:.6e}",
        })
    rows.sort(key=lambda row: float(row["p_fdr_bh"]))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    significant = [row for row in rows if float(row["p_fdr_bh"]) < alpha]
    log(f"    terms tested={len(rows)}  significant at FDR<{alpha}: {len(significant)}")
    if rows:
        best = rows[0]
        log(f"    best term: {best['GO']} ({best['name']}) FE={best['fold_enrichment']} p={best['p_uncorrected']} FDR={best['p_fdr_bh']}")
        summary.update({
            "terms_tested": len(rows),
            "significant": len(significant),
            "top_term": best["GO"],
            "top_name": best["name"],
            "top_fe": best["fold_enrichment"],
            "top_p": best["p_uncorrected"],
            "top_fdr": best["p_fdr_bh"],
        })
    return summary


def resolve_orthogroup_file(paths: PipelinePaths) -> Path:
    path = paths.root.parent / "step1" / "sources" / "orthogroups.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Orthogroup table not found: {path}. It is deposited with step 1, so "
            "run this step from inside the repository."
        )
    return path


def run_dual_level_enrichment(paths: PipelinePaths, alpha: float = ALPHA) -> Path:
    """
    Run the enrichment at both levels of analysis and write the result tables.

    The protein level is exploratory: orthologues of one family are not
    independent observations, so it counts a single ancestral gene once per
    species. The orthogroup level is the one reported in the manuscript, under
    both aggregation rules.
    """
    output_dir = paths.enrichment
    output_dir.mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    def log(message: str = "") -> None:
        print(message, flush=True)
        log_lines.append(message)

    population_file = paths.goatools / "goatools_population.txt"
    dag = GODag(str(paths.sources / "go-basic.obo"))
    associations = load_population(population_file)
    groups = load_study_files(paths.goatools)
    protein_to_orthogroup = load_protein_to_orthogroup(resolve_orthogroup_file(paths))

    log("GO enrichment - protein level and orthogroup level")
    log(f"annotated proteins (population)       : {len(associations)}")
    log(f"conservation groups found             : {len(groups)}")
    log(f"proteins mapped to an orthogroup      : {len(protein_to_orthogroup)}")

    study_accessions = [accession for group in groups.values() for accession in group]
    log(f"study accessions, distinct            : {len(set(study_accessions))}")
    if len(study_accessions) != len(set(study_accessions)):
        log(f"  WARNING: the study files hold {len(study_accessions)} lines for {len(set(study_accessions))} distinct accessions.")

    union, members = aggregate(associations, protein_to_orthogroup, "union")
    intersection, _ = aggregate(associations, protein_to_orthogroup, "intersection")
    log(f"orthogroups with >=1 annotated member : {len(union)}")

    sizes: dict[int, int] = defaultdict(int)
    for accessions in members.values():
        sizes[len(accessions)] += 1
    log("annotated members per orthogroup      : " + ", ".join(f"{size}:{count}" for size, count in sorted(sizes.items())))

    group_orthogroups = {
        name: sorted({protein_to_orthogroup[accession] for accession in accessions if accession in protein_to_orthogroup})
        for name, accessions in groups.items()
    }

    core_orthogroups = sorted(og for og in group_orthogroups.get(CORE_GROUP, []) if og in union)
    complete_core = [og for og in core_orthogroups if len(members[og]) == 4]
    accessory_proteins = [a for name, accessions in groups.items() if name != CORE_GROUP for a in accessions]
    accessory_orthogroups = sorted({og for name, ogs in group_orthogroups.items() if name != CORE_GROUP for og in ogs})
    annotated_accessory = [og for og in accessory_orthogroups if og in union]

    log(f"annotated orthogroups in the core     : {len(core_orthogroups)}")
    log(f"  of which with exactly 4 annotated proteins: {len(complete_core)}")
    log(f"annotated orthogroups outside the core: {len(annotated_accessory)}")

    (output_dir / "core_annotated_orthogroups.txt").write_text(
        "".join(f"{og}\n" for og in core_orthogroups), encoding="utf-8"
    )
    (output_dir / "accessory_annotated_orthogroups.txt").write_text(
        "".join(f"{og}\n" for og in annotated_accessory), encoding="utf-8"
    )

    summary_rows: list[dict[str, object]] = []

    def run(level: str, aggregation: str, group: str, study_ids, associations_for_level) -> None:
        log(f"-- {level} / {aggregation or '-'} / {group}")
        name = f"{level}_level__{group}.tsv" if not aggregation else f"orthogroup_{aggregation}__{group}.tsv"
        result = run_enrichment(study_ids, associations_for_level, dag, output_dir / name, log, alpha)
        summary_rows.append({"level": level, "aggregation": aggregation, "group": group, **result})

    if accessory_proteins:
        log("")
        log("### headline comparison: accessory (non-core) versus all annotated")
        run("protein", "", "accessory", accessory_proteins, associations)
        run("orthogroup", "union", "accessory", accessory_orthogroups, union)
        run("orthogroup", "intersection", "accessory", accessory_orthogroups, intersection)

    log("")
    log("### per conservation group")
    for name in sorted(groups):
        run("protein", "", name, groups[name], associations)
        run("orthogroup", "union", name, group_orthogroups[name], union)
        run("orthogroup", "intersection", name, group_orthogroups[name], intersection)

    with (output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    log("")
    log(f"results written to {output_dir}")
    (output_dir / "run_log.txt").write_text("".join(f"{line}\n" for line in log_lines), encoding="utf-8")
    return output_dir
