# 02 — Step 2: dual-level GO enrichment, wired into the default run

Covers TASK 1 of the brief. Depends on spec 01 (the study files must be correct
first, otherwise the run log's line-count warning fires and the audit is noise).

## Problem

`step2/src/opxv_analysis/enrichment.py` (45 lines) has three defects, two of
them fatal to the deposit's credibility:

1. **Wrong unit of analysis.** It runs the test with the protein as the
   observation; the manuscript reports the orthogroup. Orthologues of one family
   are not independent observations, so the protein-level test counts a single
   ancestral gene once per species. Run as deposited it yields 54 significant
   terms where the paper reports none.
2. **The study set is not a subset of the population.** `enrichment.py:30` passes
   the raw study list to `run_study`, while the population is
   `geneid2gos.keys()` — the 483 annotated proteins. For the accessory
   comparison that is 171 study ids against a population of 483, of which only 69
   are present. The extra ids inflate the study size with observations that carry
   no annotation.
3. **It does not run at all on a current environment.** GOATOOLS imports
   `multipletests` from `statsmodels.sandbox.stats.multicomp`. Recent statsmodels
   keeps the module but no longer exports the name, so `methods=["fdr_bh"]`
   raises `ImportError`. A plain `try: import …` does not detect this — the
   module imports fine.

Compounding all three: `enrichment` is a manual subcommand taking four
hand-typed paths (`run_analysis.py:126-130`), so it is not part of
`run_pipeline()` and nothing exercises it. That is how defect 2 survived.

## Do this

### 2.1 Replace `step2/src/opxv_analysis/enrichment.py`

Port the logic of the supplied `go_enrichment_dual_level.py` (repository root)
into the module. It runs the same test at two levels — protein (exploratory) and
orthogroup (the manuscript's), the latter under two aggregation rules:

- `union`: the orthogroup carries the union of its annotated members' GO terms;
- `intersection`: only the terms shared by **all** annotated members.

Unannotated members are ignored entirely, so `intersection` is across annotated
members only. Keep the functions as small units rather than one procedure:
`load_population`, `load_study_files`, `load_protein_to_orthogroup`,
`aggregate`, `run_enrichment`.

Three behaviours are load-bearing. Carry each one over, and keep the short
comment that explains why, because none of them is obvious from the code:

- **Study ∩ population before the test.** `run_enrichment` computes
  `study = sorted(set(study_ids) & population)` and logs how many ids were
  dropped. Without it, defect 2 returns.
- **A pristine deep copy of the associations into GOATOOLS.**
  `propagate_counts=True` mutates the association dict in place, so every
  subsequent run in the same process would see ancestor-closed sets. The script
  passes `copy.deepcopy({k: set(v) for k, v in assoc.items()})`.
- **The statsmodels shim**, installed *before* `from goatools… import`. Check the
  attribute, not the import:

  ```python
  if not hasattr(_sb, "multipletests"):
      _sb.multipletests = _mt.multipletests
  ```

  Creating the module when absent and patching it when present-but-empty are
  both required.

Adapt to the module's conventions while porting: `from __future__ import
annotations`, paths resolved through `PipelinePaths` (`config.default_paths()`)
rather than a `--root` argument, and no `argparse`/`main()` inside the module —
the CLI lives in `run_analysis.py`. Drop the script's `--no-self-check` block
from the module; the reference values live in this specs directory and in the
acceptance test, not in production code.

Keep the fold-enrichment computation as the script has it —
`(study_count/study_n) / (pop_count/pop_n)`, guarded against
`ZeroDivisionError` — and keep writing every tested term, not only the
significant ones. The deposited behaviour of writing a file only when something
passes FDR is itself a problem: a null result then looks like a missing file.

### 2.2 Wire it into the default flow

In `step2/run_analysis.py`:

- call the enrichment from `run_pipeline()`, after `build_goatools_inputs` (which
  produces its inputs in the same run);
- delete the four positional arguments of the `enrichment` subparser. Keep the
  subcommand itself for running the enrichment alone, but with no arguments — all
  inputs are fixed:

  | Input | Path |
  |---|---|
  | population | `step2/outputs/goatools/goatools_population.txt` |
  | study sets | `step2/outputs/goatools/goatools_study_*.txt` |
  | ontology | `step2/sources/go-basic.obo` |
  | orthogroup map | `step1/sources/orthogroups.tsv` |

  The orthogroup map is the only cross-step input; resolve it from
  `paths.root.parent / "step1" / "sources" / "orthogroups.tsv"` and fail with a
  message naming step 1 if it is missing, the way
  `step4/src/viral_filtering/services.py:19` does.
- add `enrichment` to `PipelinePaths` as `outputs / "enrichment_dual_level"` in
  `config.default_paths()`, and to whatever `setup_paths` creates, so the output
  directory is handled like every other.

The outgroup columns to drop (`Humanalphaherpesvirus1_proteins`,
`swinepox_proteins`) are already named in `step1`'s documentation and in
`config.SPECIES_MAP`'s complement; define them once in
`step2/src/opxv_analysis/config.py` rather than inline.

### 2.3 Outputs

Into `step2/outputs/enrichment_dual_level/`:

- `protein_level__<group>.tsv`, `orthogroup_union__<group>.tsv`,
  `orthogroup_intersection__<group>.tsv` — one per level, aggregation and group,
  including the `accessory` headline comparison;
- `summary.tsv` — one row per (level, aggregation, group);
- `run_log.txt`;
- `core_annotated_orthogroups.txt` (96 ids) and
  `accessory_annotated_orthogroups.txt` (29 ids). The script already computes
  both partitions but does not persist them; write them, one id per line, sorted.

## Acceptance

```bash
cd step2 && python run_analysis.py pipeline
```

Completes without error, and:

1. `outputs/enrichment_dual_level/summary.tsv` reproduces the three `GO:0005515`
   rows:

   | level / aggregation | FE | p | FDR |
   |---|---|---|---|
   | protein | 2.402 | 2.358e-09 | 2.004e-06 |
   | orthogroup, union | 2.235 | 2.007e-04 | 0.1706 |
   | orthogroup, intersection | 2.321 | 1.122e-04 | 0.0646 |

2. `run_log.txt` prints `483` annotated proteins, `125` orthogroups with
   annotation, `96` in the core, `87` of those with exactly four annotated
   proteins, and `29` outside the core.
3. Terms tested: `850` under union, `744` under intersection.
4. **No term is significant at FDR 0.05 at the orthogroup level** — this is the
   manuscript's reported result and the whole point of the task. The protein
   level will report significant terms; that file is labelled exploratory.
5. `core_annotated_orthogroups.txt` has 96 lines,
   `accessory_annotated_orthogroups.txt` has 29, and the two sets are disjoint.
6. `wc -l outputs/goatools/goatools_population.txt` still reports 483, and the
   consolidated annotation table is byte-identical to before this spec.

## Risks

- **Re-running `run_pipeline()` rewrites `step3_final_consolidated_annotations.tsv`
  and `enhanced_summary_table.tsv`.** They should come out identical, since
  nothing in their path changes. Diff them against the baseline from spec 00 and
  stop if they differ.
- The enrichment adds a GOATOOLS `GODag` load to the default run (plus one per
  level). `load_go_dag` already parses the same OBO for the annotation path;
  reuse the loaded DAG rather than parsing it again if it is cheap to thread
  through, but do not restructure the pipeline for it.
- `go_enrichment_dual_level.py` at the repository root becomes dead code once
  ported. Spec 07 removes it.
