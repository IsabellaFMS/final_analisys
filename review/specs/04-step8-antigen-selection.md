# 04 — Step 8: final antigen selection (new step)

Covers TASK 6 of the brief. Independent of specs 01–03.

This is the step that produces the manuscript's main conclusion — the two
prioritised antigens — and no code implements it. It is the one a reviewer is
most likely to ask for.

## Target result

The step must select exactly two antigens:

| Accession | Protein | Species | Orthogroup |
|---|---|---|---|
| `YP_010377182.1` | chemokine-binding protein | MPXV | `OG0000008` |
| `NP_619793.1` | CrmB (CPXV005) | CPXV | `OG0000003` |

and rule 5 must fire for `OG0000003`.

## Layout

Follow the conventions of the existing steps, closest in shape to step 6
(`run_analysis.py` + `src/<package>/` with `config`, `loaders`, `models`,
`analyzer`, `reporting`, `writer`):

```
step8/
  run_analysis.py
  README.md
  src/antigen_selection/
    __init__.py  config.py  loaders.py  models.py  analyzer.py  reporting.py  writer.py
  sources/
  outputs/
    tables/   reports/
```

`run_analysis.py` parses arguments, resolves paths and prints a summary; the
package holds the logic. Keep the module boundaries the other steps use: no
analysis in `run_analysis.py`, no file I/O in `analyzer.py`.

## Inputs

All of them already exist in the deposit. Resolve each from its sibling step,
with a fallback copy under `step8/sources/` and a `--…` override, exactly as
step 6 does for the step 1 report (`step6/README.md`, Inputs).

| Input | Path | Used for |
|---|---|---|
| B-cell epitopes at 0.2512 and 0.2012 | derived, see below | rules 1 and 5 |
| Secretion/topology consensus | `step5/outputs/tables/comparative_results_signalp6_deep_final.csv` (`;`-separated) | rule 3 |
| MHC-I candidates | `step7/outputs/mhci_final.csv` | rule 4 tie-break |
| MHC-II candidates | `step7/outputs/mhcii_final.csv` | rule 4 tie-break |
| Orthogroup membership | `step1/sources/orthogroups.tsv` | rules 2, 4, 5 |
| Conservation classes | `step1/outputs/orthogroups_consolidated_analysis_classified_protein_ids.txt` | rule 2 |
| Sequence lengths | `step2/sources/all_refseq_proteins.fasta` | rule 5 |

Lengths come from the full 832-protein FASTA, not from
`step4/outputs/fasta/safe_sequences.fasta`: the orthogroup consensus length must
be computed over **all** members, including those the host-homology filter
removed, or a truncated orthologue can look full-length against a depleted
orthogroup.

### Deriving the two epitope tables

The deposited `step6/outputs/tables/epitopes.csv` is the 0.1512 run (473
epitopes across 168 proteins). Rules 1 and 5 need the 0.2512 and 0.2012 runs.

`step6/src/bcell_analysis/analyzer.py` must not be modified and the deposited
0.1512 table must not be overwritten. Both constraints are already satisfied by
step 6's own CLI, which takes `--cutoff` and `--output-dir`:

```bash
cd step6
python run_analysis.py --cutoff 0.2512 --output-dir ../step8/outputs/bcell/0.2512
python run_analysis.py --cutoff 0.2012 --output-dir ../step8/outputs/bcell/0.2012
```

Have `step8/run_analysis.py` invoke step 6 this way as a subprocess — the
precedent is `step2/run_analysis.py:57`, which shells out to `Rscript` — and skip
the invocation when the target table already exists unless `--refresh` is passed.
Two notes:

- step 6 also writes one FASTA per protein with candidates into its output
  directory, so these subdirectories will hold ~100–170 small FASTA files each.
  That is step 6's output contract, not step 8's; leave it, and say so in the
  README.
- cross-check against the brief's independently reproduced figures: 473 / 168 /
  40 epitope occurrences across 168 / 97 / 28 proteins at 0.1512 / 0.2012 /
  0.2512. If the two new runs do not reproduce 168→97 and 40→28, stop: the
  discrepancy is in the inputs, not in step 8.

## The rules

Apply in this order, counting what each rule removes:

1. **B-cell evidence.** Keep proteins with at least one epitope in the 0.2512
   table.
2. **Conservation.** Keep proteins whose orthogroup contains all four OPXV
   species. Take this from `orthogroups.tsv` (columns other than
   `Humanalphaherpesvirus1_proteins` and `swinepox_proteins`) and cross-check
   against the step 1 class `Conserved`; they must agree, and a disagreement is a
   defect in the inputs, not something to paper over.
3. **Surface accessibility.** Keep proteins in `Group 1: High Confidence
   Consensus (SP+/THMM+)` of the step 5 comparative table. Both target antigens
   are in Group 1, as are `YP_010377181.1` and `YP_233097.1`, so this rule does
   not by itself decide the outcome — but see
   [08-open-questions.md](08-open-questions.md) §3: "per the SignalP6/DeepTMHMM
   consensus" could also admit Groups 2 and 3 (single-predictor positives).
   Implement Group 1 as the default, make the accepted groups a config constant,
   and record both readings in the README.
4. **One representative per orthogroup.** Where more than one orthologue passes,
   keep the richest epitope profile: most epitopes in the 0.2512 table, ties
   broken by combined MHC-I + MHC-II allelic coverage from the step 7 tables
   (sum of distinct alleles across that protein's peptides; `Alleles` /
   `Allele Count` columns). Resolve a remaining tie by accession, so the output
   is deterministic, and log that it happened.
5. **Truncation.** Exclude orthologues annotated as truncated or partial relative
   to the orthogroup consensus length. Where the only orthologue meeting the
   0.2512 criterion is truncated, represent the orthogroup by its full-length
   orthologue and report that one's epitopes at the intermediate stringency
   (0.2012).

### Rule 5, verified against the deposit

`OG0000003` is the case the brief names, and the numbers check out:

| Member | Species column | Length |
|---|---|---|
| `NP_619793.1` | cowpox | 355 |
| `NP_619803.1` | cowpox | 202 |
| `NP_620001.1` | cowpox | 320 |
| `NP_620005.1` | cowpox | 355 |
| `YP_010377003.1` | mpox | 349 |
| `YP_010377181.1` | mpox | 349 |
| `YP_232886.1` | vaccinia WR | 122 |
| `YP_233097.1` | vaccinia WR | 122 |
| `NP_042240.1` | varv | 349 |

`YP_233097.1` is 122 aa against a 349–355 aa consensus — truncated, as the brief
says. **The brief labels it VARV; in `orthogroups.tsv` it sits in the
`vaccinia_WR_refseq_proteins` column**, and `step6/outputs/fasta/` names it
`TNF-alpha-receptor-like_protein_Vaccinia_virus`. The length claim is right and
the species label is wrong — see
[08-open-questions.md](08-open-questions.md) §1 before this reaches a README or
the manuscript.

Define "truncated" as a config constant, not a literal: flag a member whose
length is below a fraction of the orthogroup consensus (median member length).
`122 / 352 ≈ 0.35`, while the other members sit at 0.57 and above, so a threshold
anywhere in 0.4–0.5 separates them; state the chosen value in the README and in
the run report. Treat a FASTA description containing `partial` as truncated too,
independently of length.

For `OG0000003` the rule must then substitute the full-length CPXV orthologue
`NP_619793.1` and report its 0.2012 epitopes — which is how the target result
arises, and why rule 5 cannot be implemented as a plain exclusion.

## Outputs

- `outputs/tables/final_antigens.tsv` — one row per selected antigen:
  accession, protein name, species, orthogroup, consensus length, own length,
  the stringency its epitopes are reported at, epitope count, MHC-I and MHC-II
  allele coverage, step 5 group.
- `outputs/tables/selection_trace.tsv` — one row per candidate that entered rule
  1, with a verdict column per rule (`pass`, or the rule that dropped it).
  This is the audit trail a reviewer needs; without it the two-antigen result is
  an assertion.
- `outputs/tables/selected_epitopes.tsv` — the epitopes of the selected antigens,
  carrying the step 6 columns unchanged plus the cutoff they came from.
- `outputs/reports/selection_summary.txt` — counts surviving each rule, the
  orthogroups where rule 4 or rule 5 fired, and the parameter values used.

Column names in English, matching the spelling the other steps use
(`Protein_ID`, `Classification (conservation)`).

## Acceptance

```bash
cd step8 && python run_analysis.py
cut -f1 outputs/tables/final_antigens.tsv | tail -n +2
```

1. Exactly two rows: `YP_010377182.1` and `NP_619793.1`.
2. `selection_summary.txt` records rule 5 firing for `OG0000003`, naming
   `YP_233097.1` as the excluded truncated member and `NP_619793.1` as the
   substituted full-length orthologue, with its epitopes reported at 0.2012.
3. `selection_trace.tsv` accounts for every protein that passed rule 1: the
   counts after each rule sum consistently and no protein disappears without a
   verdict.
4. Both selected accessions appear in
   `step4/outputs/reports/flagged_reintegrated_accessions.txt` (spec 03) — they
   are reincorporated duplicates, which is exactly why spec 03 matters.
5. Re-running the step produces byte-identical tables.
6. No file under `step5/outputs/`, `step6/outputs/` or `step7/outputs/` is
   modified. Verify with `git status` before committing.

## If the two antigens do not come out

Do not tune thresholds until they do. Record what the step produced, which rule
diverged, and under which reading of rules 3 and 4, and resolve it against the
manuscript. A selection pipeline reverse-fitted to its expected output is worth
nothing at review.
