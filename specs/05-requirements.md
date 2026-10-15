# 05 — requirements.txt

Covers TASK 5 of the brief. Run after spec 02, so the pins match what the fixed
code was actually verified against.

## Problem

`requirements.txt` as deposited:

```
# Steps 1, 3, 4 and 6 need nothing beyond the standard library.

# Step 2, enrichment subcommands only.
goatools>=1.4

# Steps 5 and 7.
pandas>=2.0

# Test runner.
pytest>=7.4
```

Three issues:

1. Nothing is pinned. The manuscript cites GOATOOLS v1.6.5; `>=1.4` resolves to
   whatever is current, and the enrichment results depend on it.
2. `statsmodels` is not listed at all, although GOATOOLS needs it for
   `methods=["fdr_bh"]` and the shim in spec 02 exists because of its layout.
   Resolving dependencies afresh gives no indication that the statsmodels version
   matters.
3. The comments are now wrong: after spec 02 the enrichment runs in the default
   flow, not in "enrichment subcommands only", and step 2 also needs `pandas`.

## Do this

Pin exact versions (`==`) for every direct dependency, using the versions the
fixed pipeline was verified against — take them from the baseline recorded in
spec 00 and re-check after spec 02 passes:

- `goatools==1.6.5`, the version the manuscript cites.
- `statsmodels==<verified>`, added explicitly.
- `pandas==<verified>`.
- `pytest==<verified>`, and only if tests are actually deposited: there is no
  test directory in the repository, so either add the tests or drop the
  dependency. Listing a test runner with nothing to run is one of the README/code
  mismatches this round of fixes is meant to eliminate.

Keep the per-step comments and correct them: steps 1, 3 and 4 use only the
standard library; step 2 needs `goatools`, `statsmodels` and `pandas`; steps 5
and 7 need `pandas`; step 6 needs only the standard library; step 8 (spec 04)
needs `pandas` if it uses it — say which.

Add one comment recording the statsmodels incompatibility, kept factual and
short:

```
# goatools imports multipletests from statsmodels.sandbox.stats.multicomp.
# Recent statsmodels keeps that module but no longer exports the name, so
# step 2 installs the attribute before importing goatools. Changing the
# statsmodels pin means re-checking that shim.
```

Also state the Python version the pins were verified on. The READMEs claim
"Python 3.11 or newer"; confirm it against the interpreter recorded in spec 00
and make the two agree.

## Acceptance

```bash
python -m venv /tmp/opxv && /tmp/opxv/bin/pip install -r requirements.txt
cd step2 && python run_analysis.py pipeline
```

- A fresh environment built from the pins runs spec 02's acceptance test and
  reproduces the three `GO:0005515` rows.
- `pip freeze` lists every pinned version exactly as written.
- No dependency is listed that nothing in the repository imports; check with
  `grep -rn "^import \|^from " --include=*.py .` over the step packages.
