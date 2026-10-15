# Specs — deposit fixes

Implementation specs derived from `../implementation_brief.md`. One spec per unit of
work, each with its own preconditions, exact edits, acceptance test and rollback.
They are written to be executed in order; where an order is not required the spec
says so.

The source brief states the goal: the deposit must reproduce the GO-enrichment
result reported in the manuscript, and no README statement may be contradicted by
the repository's own files.

## Order

| # | Spec | Touches | Depends on |
|---|---|---|---|
| 00 | [00-preparation.md](00-preparation.md) | nothing (snapshot, baseline) | — |
| 01 | [01-step2-study-files.md](01-step2-study-files.md) | `step2/src/opxv_analysis/analyzer.py`, `step2/outputs/goatools/` | 00 |
| 02 | [02-step2-go-enrichment.md](02-step2-go-enrichment.md) | `step2/src/opxv_analysis/enrichment.py`, `step2/run_analysis.py` | 01 |
| 03 | [03-step4-flagged-accessions.md](03-step4-flagged-accessions.md) | `step4/src/viral_filtering/` | 00 |
| 04 | [04-step8-antigen-selection.md](04-step8-antigen-selection.md) | new `step8/` | 00 |
| 05 | [05-requirements.md](05-requirements.md) | `requirements.txt` | 02 |
| 06 | [06-readme-corrections.md](06-readme-corrections.md) | `README.md`, `step2/README.md`, `step4/README.md` | 01–05 |
| 07 | [07-repository-hygiene.md](07-repository-hygiene.md) | whole repository | 01–06 |
| 08 | [08-open-questions.md](08-open-questions.md) | nothing (decisions needed) | read before 04 and 06 |

Spec 06 comes last on purpose: the READMEs must describe what the code does after
the fixes, not before.

## Reference values

Acceptance criteria shared by several specs. Every value comes from the brief and
is verifiable against the repository's own files.

| Quantity | Value |
|---|---|
| Proteins in the four OPXV proteomes | 832 |
| Proteins with at least one GO term (enrichment population) | 483 |
| Rows in `step2/outputs/tables/step3_final_consolidated_annotations.tsv` | 796 |
| Proteins with GO, by namespace | 333 BP, 230 CC, 323 MF |
| Orthogroups with at least one annotated member | 125 (96 core, 29 accessory) |
| Core orthogroups with exactly 4 annotated proteins | 87 |
| GO terms tested | 850 (union), 744 (intersection) |
| `GO:0005515`, protein level | FE 2.402, p 2.358e-09, FDR 2.004e-06 |
| `GO:0005515`, orthogroup union | FE 2.235, p 2.007e-04, FDR 0.1706 |
| `GO:0005515`, orthogroup intersection | FE 2.321, p 1.122e-04, FDR 0.0646 |
| Rows in `step4/sources/viral_hits.tsv` | 19,709 (1,078 at 12/12, 13,995 at 9/12) |
| Distinct viral proteins with at least one hit | 641 |
| Duplicates reincorporated by step 4 | 45, of which 33 twin-excluded |
| B-cell epitopes at cutoff 0.1512 | 473 across 168 proteins |
| Final antigens selected by step 8 | `YP_010377182.1`, `NP_619793.1` |

## Rules that apply to every spec

**Language.** The repository is in English — code, identifiers, comments,
file names, outputs and documentation. Keep it that way. Published column headers
(`Protein_ID`, `Classification (conservation)`, `GO_BP`, …) are data consumed by
the R step, by later steps and by the manuscript figures: do not rename them.

**Comments.** Match the density and tone already in the repository: a short
comment where a decision is not obvious from the code (see
`step4/src/viral_filtering/services.py:20` for the house style), and nothing
where the code already says it. No section banners, no restating a function
signature in prose, no changelog or progress notes in source files.

**No tool or authorship markers.** Source files, outputs and documentation carry
no references to editors, assistants or generation tooling, and no `TODO`/`FIXME`
left behind by the work. Disclosure of tooling used to prepare the deposit, if
any is required, belongs in the submission to the journal, not scattered through
the files.

**Scope.** The brief's "Do not change" list is binding:

- `step6/src/bcell_analysis/analyzer.py` — the peak-and-extension implementation
  is correct and independently reproduced.
- the MHC-II three-cutoff consolidation in step 7 — redundant, not wrong.
- any numerical result under `step2/outputs/`, `step5/outputs/`,
  `step6/outputs/` or `step7/outputs/` other than the files a spec regenerates.

**Verification before assertion.** No spec is complete until its acceptance test
has actually been run and its output recorded. A number that cannot be
reproduced from the deposit does not go into a README; it goes into
[08-open-questions.md](08-open-questions.md) for the author to resolve.
