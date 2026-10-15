# Deposit fixes — Orthopoxvirus comparative proteomics pipeline

Implementation brief. Each task has an exact location, the exact current text where
applicable, and an acceptance test. Do not change anything outside these tasks.

**Repository**: the seven-step deposit (`step1/` … `step7/`, plus root `README.md`
and `requirements.txt`). Associated manuscript: *Comparative Proteomic and
Immunoinformatic Mapping of Orthopoxviruses*, submitted to PLOS Computational Biology.

**Why this matters**: the deposit currently produces a GO-enrichment result that
contradicts the one reported in the manuscript, and several README statements are
contradicted by the repository's own files. Reviewers are expected to download and
run this.

---

## Reference values

Every number below was verified against the repository's own files. Use them as
acceptance criteria.

| Quantity | Value |
|---|---|
| Proteins in the four OPXV proteomes | 832 |
| Proteins with at least one GO term (enrichment population) | 483 |
| Rows in `step2/outputs/tables/step3_final_consolidated_annotations.tsv` | 796 |
| Proteins with GO, by namespace | 333 BP, 230 CC, 323 MF |
| Orthogroups with at least one annotated member | 125 |
| …in the conserved core | 96 |
| …outside the core | 29 |
| Core orthogroups with exactly 4 annotated proteins | 87 |
| GO terms tested (union rule) | 850 |
| GO terms tested (intersection rule) | 744 |
| Protein level, `GO:0005515` | FE 2.402, p 2.358e-09, FDR 2.004e-06 |
| Orthogroup level, union, `GO:0005515` | FE 2.235, p 2.007e-04, FDR 0.1706 |
| Orthogroup level, intersection, `GO:0005515` | FE 2.321, p 1.122e-04, FDR 0.0646 |
| Rows in `step4/sources/viral_hits.tsv` | 19,709 |
| …with 12/12 identical residues | 1,078 |
| …with 9/12 | 13,995 |
| Distinct viral proteins with at least one hit | 641 |
| Duplicate accessions reincorporated by step 4 | 45 |
| …identical to an accession excluded for host homology | 33 |

---

## TASK 1 — step 2: replace the GO enrichment and wire it into the pipeline

**Problem.** `step2/src/opxv_analysis/enrichment.py` runs the test with the
**protein** as the unit of analysis. The manuscript reports the test with the
**orthogroup** as the unit. Running the deposit as it stands yields 54 significant
terms where the paper reports none.

Two additional defects in the same path:

1. The study sets are not subsets of the population. The study files hold all 832
   proteins; the population holds only the 483 that carry a GO term. For the
   accessory (non-core) comparison that is 171 distinct study ids against a
   population of 483, of which only 69 are in the population.
2. `enrichment` is a manual subcommand with four hand-typed paths
   (`run_analysis.py enrichment <study> <population> <obo> <output>`), which is how
   defect 1 went unnoticed.

**Do this.**

- Replace the contents of `step2/src/opxv_analysis/enrichment.py` with the logic of
  the supplied `go_enrichment_dual_level.py` (delivered alongside this brief). It
  runs the test at both levels — protein (exploratory) and orthogroup (the
  manuscript's), the latter under two aggregation rules, union and intersection.
- Keep its three behaviours, they are load-bearing:
  - the study set is intersected with the population before the test;
  - a pristine deep copy of the associations is passed to GOATOOLS, because
    `propagate_counts=True` mutates the association dict in place;
  - the statsmodels shim. GOATOOLS imports `multipletests` from
    `statsmodels.sandbox.stats.multicomp`; recent statsmodels keeps that module but
    no longer exports the name, so a plain `try: import …` does **not** detect it.
    Check for the attribute, not the import. **Without this the current code raises
    `ImportError` and the step does not run at all.**
- Wire it into the default flow of `step2/run_analysis.py`. No hand-typed paths: the
  inputs are always `step2/outputs/goatools/goatools_population.txt`,
  `step2/outputs/goatools/goatools_study_*.txt`, `step2/sources/go-basic.obo` and
  `step1/sources/orthogroups.tsv`. Keep a subcommand if you like, but the full run
  must produce the enrichment without arguments.
- Write the outputs to `step2/outputs/enrichment_dual_level/`: one TSV per level,
  aggregation and group, plus `summary.tsv` and `run_log.txt`.
- Additionally write two partition lists (the script already computes them):
  `step2/outputs/enrichment_dual_level/core_annotated_orthogroups.txt` (96 ids) and
  `accessory_annotated_orthogroups.txt` (29 ids).

**Acceptance.** `python run_analysis.py` in `step2/` completes without error, and
`summary.tsv` reproduces the three `GO:0005515` rows in the reference table above.
The run log prints `483`, `125`, `96`, `87` and `29`.

---

## TASK 2 — step 2: fix the study-file writer

**Problem.** The partially-conserved study files repeat each accession once per
species in the combination. Observed:

| File | Lines | Distinct |
|---|---|---|
| `goatools_study_CPXV_plus_VACV_plus_VARV.txt` | 150 | 50 |
| `goatools_study_CPXV_plus_MPXV_plus_VACV.txt` | 90 | 30 |
| `goatools_study_CPXV_plus_MPXV_plus_VARV.txt` | 21 | 7 |
| `goatools_study_CPXV_plus_VACV.txt` | 16 | 8 |
| `goatools_study_MPXV_plus_VACV.txt` | 12 | 6 |
| `goatools_study_VACV_plus_VARV.txt` | 8 | 4 |
| `goatools_study_CPXV_plus_MPXV.txt` | 4 | 2 |

The multiplier is the number of species in the combination. `goatools_study_core.txt`
and the species-restricted files are not affected. There is **no overlap between
files** — the conservation classes are mutually exclusive, which is correct.

GOATOOLS deduplicates internally, so the test result is unaffected, but the files are
wrong and they invalidate any audit by line count.

**Do this.** Find the writer that emits these files (and the matching
`step2/sources/list_*.txt`) and make it emit one line per protein. Then regenerate
them.

**Acceptance.** Every `goatools_study_*.txt` has as many lines as distinct
accessions. The eleven files together hold 832 distinct accessions, none of them from
the `Humanalphaherpesvirus1_proteins` or `swinepox_proteins` columns of
`step1/sources/orthogroups.tsv`.

---

## TASK 3 — step 4: persist the flagged accessions

**Problem.** `step4/README.md` states that 33 of the 45 reincorporated duplicates are
identical to accessions excluded for host homology, and that they are flagged. The
code computes exactly that — `analyzer.py` builds
`DuplicateVerdict(..., twin_excluded=kept_instead in exclusion_ids, ...)` and
`reporting.py` prints `of which identical to an excluded accession: 33` — but the
list itself is never written, so nothing is actually flagged.

The manuscript's Results section states that these accessions are flagged in the
deposited data. Without the file that sentence has no support.

**Do this.** Write `step4/outputs/reports/flagged_reintegrated_accessions.txt` (TSV
with a header) from the existing `ReincorporationSummary.verdicts`, one row per
verdict with `twin_excluded` true. Columns: `protein_id`, `kept_instead`, `pathogen`,
`source_index`.

**Acceptance.** The file has 33 data rows. `YP_010377182.1` and `NP_619793.1` are
among them — these are the two antigens the manuscript prioritizes, and a reviewer
will look for exactly those.

---

## TASK 4 — README corrections

These are factual errors: the repository's own files contradict them.

### `README.md` (root)

| Line | Current | Replace with |
|---|---|---|
| 3 | `Seven analysis steps used to select conserved…` | Leave as is. The *manuscript* says "six sequential modules" and will be corrected instead. Do not change the README here. |
| 29 | `The proteome of each virus is the background population.` | `The background population is the set of 483 proteins carrying at least one GO term, pooled across the four proteomes.` The code uses `geneid2gos.keys()`, i.e. the contents of `goatools_population.txt`. |
| 41 | `…share identical 12-residue stretches with the human proteome…` | `…contain a 12-residue window with 9 or more residues identical to a human 12-mer…` Only 1,078 of 19,709 hit rows are exact matches; 13,995 are 9/12. Also record the parameters of the external k-mer run, since the Maguire et al. scripts are not redistributed and this threshold is the only thing that makes the filter reproducible. |
| 55 | `…extracts contiguous regions above the score cutoff, keeping those between 20 and 30 residues.` | `…masks runs of more than three consecutive sub-threshold residues, anchors each candidate at the highest remaining score and grows it symmetrically, keeping regions of 20 to 30 residues.` `step6/README.md` already describes this correctly — align the root summary to it, not the other way round. |
| 68 | data-flow row: `orthogroups.tsv, orthogroups_unassigned_genes.tsv \| OrthoFinder` | The deposited `orthogroups.tsv` has 339 rows (`OG0000000`–`OG0000338`, 1,059 proteins); OrthoFinder's file of that name has 211. Either rename it to `orthogroups_merged.tsv` and update every reference, or state in this row that it is the concatenation. `step1/README.md` already explains the 128-id overlap; the root table does not. |
| 112–122 | versions `not recorded` | OrthoFinder v2.5.5 (MCL inflation 1.7); GOATOOLS v1.6.5. Leave DeepTMHMM as not recorded if it genuinely was not, but say so explicitly. InterProScan reads `v5.75-106.019` here and `v5.75-106.0` in the manuscript — determine which is right and make both agree. |

Also add to the root README:

- The 36 proteins (37–86 aa) absent from
  `step2/outputs/tables/step3_final_consolidated_annotations.tsv`: neither
  InterProScan nor eggNOG-mapper returned anything for them. The table therefore has
  796 rows, not 832, while the manuscript's denominator is 832, counting them as
  unannotated.
- That pruning redundant parent GO terms upstream is statistically neutral given
  `propagate_counts=True`. Verified: intersecting the pruned sets and intersecting
  ancestor-closed sets both give 744 terms tested and FDR 0.0646.
- That the run log reports 41 obsolete GO ids absent from the deposited
  `go-basic.obo` (release 2024-01-17). GOATOOLS skips them without error.

### `step2/README.md`

| Line | Current | Replace with |
|---|---|---|
| 1 | `# Step 2: functional annotation and GO enrichment` | `# Step 2: functional annotation and GO term distribution` — the reported result is null, so "enrichment" in the title misleads. |
| 11 | `…with the proteome of each virus as the background population.` | Same replacement as root line 29, worded identically. |
| 57 | `Of the 796 proteins, 241 carry at least one GO term: 150 with a biological process, 168 with a cellular component and 59 with a molecular function.` | `Of the 796 proteins, 483 carry at least one GO term: 333 with a biological process, 230 with a cellular component and 323 with a molecular function.` These are stale numbers from an earlier run; the table in the same directory gives the correct ones. |

Also add a description of the two-level analysis and the two aggregation rules.

### `step4/README.md`

- Same `identical 12-residue` correction as root line 41.
- Record the external k-mer run's parameters.

### `step6/README.md`

- **No changes.** It is correct and is the reference for the root summary.

---

## TASK 5 — `requirements.txt`

- Pin exact versions. The manuscript cites GOATOOLS v1.6.5; the file currently says
  `goatools>=1.4`.
- Add a comment about the statsmodels incompatibility described in Task 1, so anyone
  resolving dependencies afresh knows why the shim exists.

---

## TASK 6 — step 8: final antigen selection (new)

**Problem.** The manuscript's *Final Antigen Selection* section describes precise
rules, and no code implements them. This is the step that produces the paper's main
conclusion and the one most likely to be requested at review.

**Do this.** Create `step8/` following the conventions of the other steps
(`run_analysis.py`, `src/`, `sources/`, `outputs/`, `README.md`). Starting from the
step 5, 6 and 7 tables, apply in order:

1. keep proteins with at least one B-cell epitope at BepiPred cutoff 0.2512;
2. restrict to proteins conserved across all four OPXV species;
3. require surface accessibility per the SignalP6 / DeepTMHMM consensus;
4. where more than one orthologue of the same orthogroup passes, pick the one with
   the richest epitope profile — more conserved B-cell epitopes, ties broken by
   combined MHC class I and class II allelic coverage;
5. exclude orthologues annotated as truncated or partial relative to the orthogroup
   consensus length; where the only orthologue meeting the high-stringency B-cell
   criterion is truncated, represent the orthogroup by its full-length orthologue and
   report its epitopes at the intermediate stringency (0.2012).

**Acceptance.** The step selects exactly two antigens: `YP_010377182.1` (MPXV
chemokine-binding protein) and `NP_619793.1` (CPXV CrmB). Rule 5 must fire for the
CrmB orthogroup `OG0000003`, whose only member with an epitope at 0.2512 is
`YP_233097.1` (VARV, 122 aa against a 349–355 aa consensus).

---

## Do not change

- `step6/src/bcell_analysis/analyzer.py`. The peak-and-extension implementation is
  correct and matches the manuscript; an independent reimplementation reproduced its
  473 / 168 / 40 occurrences across 168 / 97 / 28 proteins at the three cutoffs.
- The MHC-II three-cutoff consolidation in step 7. Filtering at 5, 7.5 and 10 and
  concatenating with `deduplicate=True` yields the same set as filtering at 10 alone,
  because the first two are subsets of the third. It is redundant, not wrong.
- Any numerical result in `step2/outputs/`, `step5/outputs/`, `step6/outputs/` or
  `step7/outputs/` other than the files these tasks regenerate.
