#!/usr/bin/env python3
"""
GO enrichment for the Orthopoxvirus comparative proteomics pipeline — step 2.

Runs the SAME enrichment test at TWO levels of analysis:

  (a) PROTEIN level  - one observation per protein. Exploratory only: orthologues
                       of the same family are not independent observations, so
                       this test counts a single ancestral gene once per species.
  (b) ORTHOGROUP level - one observation per orthogroup. This is the test reported
                       in the manuscript. Two aggregation rules are applied:
                         union        : the orthogroup carries the union of the GO
                                        terms of its annotated members
                         intersection : the orthogroup carries only the GO terms
                                        shared by ALL its annotated members

Statistics in every case: two-sided Fisher's exact test (GOATOOLS default
pvalcalc='fisher_scipy_stats') with Benjamini-Hochberg FDR correction and counts
propagated to ancestor terms.

Inputs (all already present in the deposit):
  step1/sources/orthogroups.tsv                      protein -> orthogroup map
  step2/sources/go-basic.obo                         ontology
  step2/outputs/goatools/goatools_population.txt     id2gos, one line per annotated protein
  step2/outputs/goatools/goatools_study_*.txt        one accession per line, per conservation group

Outputs (written to step2/outputs/enrichment_dual_level/):
  protein_level__<group>.tsv
  orthogroup_union__<group>.tsv
  orthogroup_intersection__<group>.tsv
  summary.tsv                 one row per (level, aggregation, group)
  run_log.txt

Usage:
  python go_enrichment_dual_level.py --root /path/to/deposit [--alpha 0.05] [--no-self-check]

Requires: goatools>=1.6.5
"""

from __future__ import annotations

import argparse
import copy
import csv
import os
import sys
import types
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# statsmodels compatibility shim.
# GOATOOLS imports multipletests from statsmodels.sandbox.stats.multicomp,
# a path removed in recent statsmodels releases. Recreate it before importing
# GOATOOLS so that methods=["fdr_bh"] keeps working.
# ---------------------------------------------------------------------------
try:
    import statsmodels.stats.multitest as _mt

    try:
        import statsmodels.sandbox.stats.multicomp as _sb  # may exist but be empty
    except Exception:
        _sb = types.ModuleType("statsmodels.sandbox.stats.multicomp")
        sys.modules["statsmodels.sandbox.stats.multicomp"] = _sb
    if not hasattr(_sb, "multipletests"):
        _sb.multipletests = _mt.multipletests
except Exception:  # pragma: no cover - statsmodels absent
    pass

from goatools.go_enrichment import GOEnrichmentStudy  # noqa: E402
from goatools.obo_parser import GODag  # noqa: E402


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_population(path: Path) -> dict[str, set[str]]:
    """goatools_population.txt: '<accession>\\t<GO:x;GO:y;...>' per line."""
    assoc: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            acc = parts[0].strip()
            gos = {g.strip() for g in parts[1].split(";") if g.strip()}
            if acc and gos:
                assoc[acc] = gos
    return assoc


def load_study_files(d: Path) -> dict[str, list[str]]:
    """One file per conservation group; one accession per line."""
    groups: dict[str, list[str]] = {}
    for p in sorted(d.glob("goatools_study_*.txt")):
        name = p.stem[len("goatools_study_"):]
        with p.open(encoding="utf-8") as fh:
            accs = [l.strip() for l in fh if l.strip()]
        groups[name] = accs
    return groups


def load_protein_to_orthogroup(path: Path, drop_columns: tuple[str, ...]) -> dict[str, str]:
    """
    orthogroups.tsv as deposited: header 'Orthogroup' + one column per proteome.
    Outgroup columns are dropped so the mapping covers only the four OPXV species,
    matching the conservation classification produced by step 1.
    """
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        keep = [i for i, h in enumerate(header)
                if i > 0 and h.strip() not in drop_columns]
        for row in rdr:
            if not row:
                continue
            og = row[0].strip()
            for i in keep:
                if i >= len(row):
                    continue
                for acc in row[i].split(","):
                    acc = acc.strip()
                    if acc:
                        mapping[acc] = og
    return mapping


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def aggregate(assoc: dict[str, set[str]],
              prot2og: dict[str, str],
              rule: str) -> tuple[dict[str, set[str]], dict[str, list[str]]]:
    """
    Collapse a protein-level id2gos into an orthogroup-level id2gos.
    Only annotated proteins contribute. Unannotated members are ignored entirely,
    so 'intersection' is the intersection across ANNOTATED members only.
    """
    if rule not in {"union", "intersection"}:
        raise ValueError(f"unknown aggregation rule: {rule}")

    members: dict[str, list[str]] = defaultdict(list)
    for acc in assoc:
        og = prot2og.get(acc)
        if og is None:
            continue
        members[og].append(acc)

    og_assoc: dict[str, set[str]] = {}
    for og, accs in members.items():
        sets = [assoc[a] for a in accs]
        if rule == "union":
            terms: set[str] = set().union(*sets)
        else:
            terms = set(sets[0]).intersection(*sets[1:]) if len(sets) > 1 else set(sets[0])
        if terms:
            og_assoc[og] = terms
    return og_assoc, members


# ---------------------------------------------------------------------------
# Enrichment
# ---------------------------------------------------------------------------
def run_enrichment(study_ids: list[str],
                   assoc: dict[str, set[str]],
                   dag: GODag,
                   alpha: float,
                   out_path: Path,
                   log) -> dict:
    """
    One enrichment run. The study set is ALWAYS restricted to the population
    before the test: passing ids absent from the population inflates the study
    size with observations that carry no annotation and distorts the p-values.

    GOATOOLS mutates the association dict in place when propagate_counts=True,
    so a pristine copy is kept for any downstream bookkeeping.
    """
    population = set(assoc)
    study = sorted(set(study_ids) & population)
    dropped = sorted(set(study_ids) - population)

    log(f"    population={len(population)}  study_in={len(study_ids)}  "
        f"study_used={len(study)}  dropped_not_in_population={len(dropped)}")

    if not study:
        log("    study is empty after restriction to the population; skipped")
        return {"population": len(population), "study_in": len(study_ids),
                "study_used": 0, "dropped": len(dropped),
                "terms_tested": 0, "significant": 0,
                "top_term": "", "top_name": "", "top_fe": "",
                "top_p": "", "top_fdr": ""}

    working = copy.deepcopy({k: set(v) for k, v in assoc.items()})

    goea = GOEnrichmentStudy(
        list(working.keys()),
        working,
        dag,
        propagate_counts=True,
        alpha=alpha,
        methods=["fdr_bh"],
    )
    results = goea.run_study(study)

    rows = []
    for r in results:
        ratio_s = getattr(r, "ratio_in_study", (0, 0))
        ratio_p = getattr(r, "ratio_in_pop", (0, 0))
        try:
            fe = (ratio_s[0] / ratio_s[1]) / (ratio_p[0] / ratio_p[1])
        except ZeroDivisionError:
            fe = float("nan")
        rows.append({
            "GO": r.GO,
            "name": r.name,
            "namespace": getattr(r, "NS", ""),
            "enrichment": getattr(r, "enrichment", ""),
            "study_count": ratio_s[0], "study_n": ratio_s[1],
            "pop_count": ratio_p[0], "pop_n": ratio_p[1],
            "fold_enrichment": f"{fe:.6g}",
            "p_uncorrected": f"{r.p_uncorrected:.6e}",
            "p_fdr_bh": f"{r.p_fdr_bh:.6e}",
        })

    rows.sort(key=lambda d: float(d["p_fdr_bh"]))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    sig = [r for r in rows if float(r["p_fdr_bh"]) < alpha]
    log(f"    terms tested={len(rows)}  significant at FDR<{alpha}: {len(sig)}")
    if rows:
        t = rows[0]
        log(f"    best term: {t['GO']} ({t['name']}) "
            f"FE={t['fold_enrichment']} p={t['p_uncorrected']} FDR={t['p_fdr_bh']}")

    top = rows[0] if rows else {}
    return {"population": len(population), "study_in": len(study_ids),
            "study_used": len(study), "dropped": len(dropped),
            "terms_tested": len(rows), "significant": len(sig),
            "top_term": top.get("GO", ""), "top_name": top.get("name", ""),
            "top_fe": top.get("fold_enrichment", ""),
            "top_p": top.get("p_uncorrected", ""),
            "top_fdr": top.get("p_fdr_bh", "")}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path,
                    help="deposit root (the directory containing step1/, step2/, ...)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--core-group", default="core",
                    help="name of the study file holding the conserved core "
                         "(default: core, i.e. goatools_study_core.txt)")
    ap.add_argument("--outgroup-columns", nargs="*",
                    default=["Humanalphaherpesvirus1_proteins", "swinepox_proteins"],
                    help="orthogroups.tsv columns to drop (outgroups)")
    ap.add_argument("--no-self-check", action="store_true")
    args = ap.parse_args()

    root: Path = args.root.resolve()
    out_dir = root / "step2" / "outputs" / "enrichment_dual_level"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run_log.txt"
    log_fh = log_path.open("w", encoding="utf-8")

    def log(msg: str = "") -> None:
        print(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    obo = root / "step2" / "sources" / "go-basic.obo"
    pop_file = root / "step2" / "outputs" / "goatools" / "goatools_population.txt"
    study_dir = root / "step2" / "outputs" / "goatools"
    og_file = root / "step1" / "sources" / "orthogroups.tsv"

    for p in (obo, pop_file, og_file):
        if not p.is_file():
            log(f"MISSING INPUT: {p}")
            return 2

    log("=" * 78)
    log("GO enrichment — protein level and orthogroup level")
    log("=" * 78)

    dag = GODag(str(obo))
    assoc = load_population(pop_file)
    groups = load_study_files(study_dir)
    prot2og = load_protein_to_orthogroup(og_file, tuple(args.outgroup_columns))

    log(f"annotated proteins (population)        : {len(assoc)}")
    log(f"conservation groups found             : {len(groups)}")
    log(f"proteins mapped to an orthogroup      : {len(prot2og)}")
    _all = [a for v in groups.values() for a in v]
    log(f"study accessions, distinct            : {len(set(_all))}")
    if len(_all) != len(set(_all)):
        log(f"  WARNING: the study files hold {len(_all)} lines for {len(set(_all))} "
            "distinct accessions. Some accessions are repeated, which makes any "
            "line-count audit of these files wrong. Counts below use distinct ids.")
        for _name, _v in sorted(groups.items()):
            if len(_v) != len(set(_v)):
                log(f"    {_name}: {len(_v)} lines, {len(set(_v))} distinct")

    og_union, members = aggregate(assoc, prot2og, "union")
    og_inter, _ = aggregate(assoc, prot2og, "intersection")
    log(f"orthogroups with >=1 annotated member : {len(og_union)}")

    sizes = defaultdict(int)
    for og, accs in members.items():
        sizes[len(accs)] += 1
    log("annotated members per orthogroup      : " +
        ", ".join(f"{k}:{v}" for k, v in sorted(sizes.items())))

    # map each conservation group to its orthogroups
    group_ogs: dict[str, list[str]] = {}
    for name, accs in groups.items():
        ogs = {prot2og[a] for a in accs if a in prot2og}
        group_ogs[name] = sorted(ogs)

    core = args.core_group
    if core in groups:
        _core_ogs = {og for og in group_ogs[core] if og in og_union}
        _core_four = [og for og in _core_ogs if len(members[og]) == 4]
        log(f"annotated orthogroups in the core     : {len(_core_ogs)}")
        log(f"  of which with exactly 4 annotated proteins: {len(_core_four)}")
        log("  (this is the ratio quoted in the cover letter: "
            f"{len(_core_four)} of {len(_core_ogs)})")
        accessory_prot = [a for g, accs in groups.items() if g != core for a in accs]
        accessory_og = sorted({og for g, ogs in group_ogs.items() if g != core
                               for og in ogs})
    else:
        log(f"WARNING: core group '{core}' not found; accessory set not built")
        accessory_prot, accessory_og = [], []

    # ---------------- self-check ----------------
    if not args.no_self_check:
        log()
        log("self-check against the values reported in the manuscript")
        checks = [
            ("annotated proteins", len(assoc), 483),
            ("orthogroups with annotation", len(og_union), 125),
            ("annotated orthogroups outside the core",
             len([o for o in accessory_og if o in og_union]), 29),
        ]
        for label, got, want in checks:
            flag = "OK " if got == want else "DIVERGES"
            log(f"  [{flag}] {label}: got {got}, manuscript reports {want}")
        if any(g != w for _, g, w in checks):
            log("  At least one value diverges. Inspect before using these results.")
        log()

    summary_rows = []

    def do(level: str, aggregation: str, group: str,
           study_ids: list[str], a: dict[str, set[str]]) -> None:
        log(f"-- {level} / {aggregation or '-'} / {group}")
        fn = (f"{level}_level__{group}.tsv" if not aggregation
              else f"orthogroup_{aggregation}__{group}.tsv")
        res = run_enrichment(study_ids, a, dag, args.alpha, out_dir / fn, log)
        summary_rows.append({"level": level, "aggregation": aggregation,
                             "group": group, **res})

    # headline comparison: accessory versus everything annotated
    if accessory_prot:
        log()
        log("### headline comparison: accessory (non-core) versus all annotated")
        do("protein", "", "accessory", accessory_prot, assoc)
        do("orthogroup", "union", "accessory", accessory_og, og_union)
        do("orthogroup", "intersection", "accessory", accessory_og, og_inter)

    # every conservation group, both levels
    log()
    log("### per conservation group")
    for name in sorted(groups):
        do("protein", "", name, groups[name], assoc)
        do("orthogroup", "union", name, group_ogs[name], og_union)
        do("orthogroup", "intersection", name, group_ogs[name], og_inter)

    with (out_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    log()
    log(f"results written to {out_dir}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
