# Review material — not part of the deposit

Working documents for the round of corrections applied to the deposit. They
record why the repository changed and how each change was verified, and they are
kept under version control for that reason.

**Do not upload this directory with the deposit.** It describes defects, not the
method, and none of it is a result. What gets deposited is the root `README.md`,
`requirements.txt` and `step1/` to `step8/`.

- `implementation_brief.md` — the brief the corrections were taken from, with the
  reference values each fix was accepted against.
- `specs/` — one spec per unit of work: the problem, the verified root cause, the
  change, the acceptance test and the rollback. `specs/08-open-questions.md`
  holds the five points that need an answer from the authors; they are listed in
  the summary below.

## What was done

| Spec | Change | Status |
|---|---|---|
| 00 | Version control initialised, baseline recorded, pre-fix defect reproduced | done |
| 01 | Study files write one line per protein, 832 distinct accessions | done |
| 02 | GO test runs at protein and orthogroup level, in the default run | done |
| 03 | The 33 flagged reincorporated accessions are written | done |
| 04 | New `step8/`, selecting `NP_619793.1` and `YP_010377182.1` | withdrawn |
| 05 | Dependencies pinned to the verified versions | done |
| 06 | README statements corrected against the files | done |
| 07 | Stripped-docstring residue, machine paths, duplicated doc removed | done |

Each spec's acceptance test was run; the results are in the commit messages and
in the step outputs.

Spec 04 was implemented, verified against its acceptance criteria and then
withdrawn: the authors did not want `step8/` in the deposit. The implementation
is recoverable from the commit that removed it. Nothing else depended on it, so
the only other change was undoing the references to step 8 in the root README,
`step4/README.md` and `requirements.txt`. The manuscript still describes the
selection rules, so a reviewer may ask for the code that applies them.
