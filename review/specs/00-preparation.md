# 00 — Preparation: snapshot and baseline

No file in the repository is modified by this spec. It exists because specs 01,
02 and 03 regenerate deposited outputs, and there is currently no way to undo
that.

## Why

The working directory is **not a git repository** (no `.git/`). Several specs
overwrite files under `step2/outputs/` and `step4/outputs/` that are part of the
deposit and whose current contents are referenced by the manuscript. Without a
snapshot, a regeneration that turns out wrong is unrecoverable.

## Do this

1. Initialise version control at the repository root and commit the deposit as
   it stands, so every later change is reviewable as a diff:

   ```bash
   git init
   git add -A
   git commit -m "Deposit as submitted, before the review fixes"
   ```

   `.gitignore` already exists; check that it does not exclude deposited outputs
   before committing, because those are exactly the files that must be tracked.

2. If version control is not wanted, copy the whole directory to a sibling
   backup instead, and record the date of the copy.

3. Record the baseline so later numbers can be compared against it. Run from the
   repository root:

   ```bash
   wc -l step2/outputs/tables/step3_final_consolidated_annotations.tsv \
         step2/outputs/goatools/goatools_population.txt \
         step2/outputs/goatools/goatools_study_*.txt \
         step4/sources/viral_hits.tsv \
         step6/outputs/tables/epitopes.csv
   python -c "import goatools, statsmodels, pandas; print(goatools.__version__, statsmodels.__version__, pandas.__version__)"
   python --version
   ```

   Keep the output outside the repository (it is a working note, not a deposit
   artifact).

4. Confirm the environment reproduces the defect spec 02 fixes, so the fix is
   demonstrably a fix:

   ```bash
   cd step2
   python run_analysis.py enrichment outputs/goatools/goatools_study_core.txt \
       outputs/goatools/goatools_population.txt sources/go-basic.obo /tmp/baseline.tsv
   ```

   Expected on a current statsmodels: `ImportError` on
   `statsmodels.sandbox.stats.multicomp.multipletests`. If instead it completes,
   record how many significant terms it reports — the brief observes 54 where the
   manuscript reports none.

## Acceptance

- The deposit is recoverable: either `git log` shows the baseline commit with a
  clean tree, or a dated backup copy exists.
- The baseline line counts and package versions are recorded.
- The pre-fix behaviour of the enrichment path is recorded (ImportError, or the
  number of significant terms).

## Notes

`step6/sources/0.1512_raw_output.csv`, `step2/sources/all_refseq_proteins.fasta`
and `step2/sources/go-basic.obo` are large inputs that no spec modifies. If the
backup route is taken rather than git, they can be excluded from the copy, but
then the copy is not a complete deposit — label it accordingly.
