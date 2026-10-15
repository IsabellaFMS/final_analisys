# 01 — Step 2: one line per protein in the study files

Covers TASK 2 of the brief.

## Problem

The partially-conserved study files under `step2/outputs/goatools/` repeat each
accession once per species in the combination:

| File | Lines | Distinct |
|---|---|---|
| `goatools_study_CPXV_plus_VACV_plus_VARV.txt` | 150 | 50 |
| `goatools_study_CPXV_plus_MPXV_plus_VACV.txt` | 90 | 30 |
| `goatools_study_CPXV_plus_MPXV_plus_VARV.txt` | 21 | 7 |
| `goatools_study_CPXV_plus_VACV.txt` | 16 | 8 |
| `goatools_study_MPXV_plus_VACV.txt` | 12 | 6 |
| `goatools_study_VACV_plus_VARV.txt` | 8 | 4 |
| `goatools_study_CPXV_plus_MPXV.txt` | 4 | 2 |

GOATOOLS deduplicates internally, so no published number changes. The files
themselves are wrong, and any audit by line count is invalid.

## Root cause

Verified against the files, and it differs from the brief's assumption.

The `step2/sources/list_*.txt` inputs are **already correct** — one line per
protein, 832 distinct accessions across the eleven files, no duplicates in any
file. They need no fix and must not be regenerated:

```
list_core 661   list_exclusive_CPXV 32   list_exclusive_VACV 18
list_exclusive_VARV 14   list_partial_CPXV_MPXV 2
list_partial_CPXV_MPXV_VACV 30   list_partial_CPXV_MPXV_VARV 7
list_partial_CPXV_VACV 8   list_partial_CPXV_VACV_VARV 50
list_partial_MPXV_VACV 6   list_partial_VACV_VARV 4        total 832
```

The multiplication is introduced downstream.
`load_conservation_groups` (`step2/src/opxv_analysis/loader.py:38`) deliberately
emits one row per `(protein, species)` pair, because the `Species` column is what
the group summary and the heatmap need. `build_goatools_inputs`
(`step2/src/opxv_analysis/analyzer.py:175`) then writes every row of that frame:

```python
proteins = conservation[conservation["conservation_group"] == group]["protein_id"].tolist()
```

while the sibling function `build_group_summary` in the same module already
guards against this with `drop_duplicates(subset=["protein_id"])`. Note that
`build_goatools_inputs` *does* deduplicate when deciding which groups exist
(`conservation.drop_duplicates(subset=["protein_id"])["conservation_group"]`) but
not when listing a group's members — which is why the group set is right and the
line counts are not.

## Do this

In `step2/src/opxv_analysis/analyzer.py`, inside `build_goatools_inputs`,
deduplicate by protein before writing each group file, matching the idiom
`build_group_summary` already uses. Two constraints:

- Preserve the existing order of accessions within a file (keep first
  occurrence); do not sort. Sorting would produce a large diff in files whose
  content is otherwise unchanged.
- End each file with a trailing newline. The current writer uses
  `"\n".join(proteins)`, which leaves the last line unterminated; `wc -l`
  therefore under-reports by one. Switch to one `f"{protein}\n"` per line, as the
  population writer two lines above already does.

Do not touch `load_conservation_groups`: the per-species expansion is required by
`build_group_summary` and by the heatmap, and changing it would alter published
outputs that no spec authorises changing.

Regenerate afterwards:

```bash
cd step2 && python run_analysis.py pipeline
```

`pipeline` rather than the default `full`, so the R heatmap step is not required
to be installed for this spec. The heatmap reads the annotation tables, which
this spec does not change.

## Acceptance

```bash
cd step2/outputs/goatools
for f in goatools_study_*.txt; do
  printf '%s %s %s\n' "$f" "$(grep -c . "$f")" "$(sort -u "$f" | grep -c .)"
done
cat goatools_study_*.txt | sort -u | grep -c .
```

- Every file reports the same count twice.
- The eleven files together hold **832** distinct accessions.
- No accession comes from the `Humanalphaherpesvirus1_proteins` or
  `swinepox_proteins` columns of `step1/sources/orthogroups.tsv`:

  ```bash
  python - <<'PY'
  import csv, pathlib
  root = pathlib.Path("..") / ".." / ".."
  rows = list(csv.reader((root / "step1/sources/orthogroups.tsv").open(encoding="utf-8"), delimiter="\t"))
  header = rows[0]
  drop = {"Humanalphaherpesvirus1_proteins", "swinepox_proteins"}
  idx = [i for i, h in enumerate(header) if h.strip() in drop]
  outgroup = {a.strip() for r in rows[1:] for i in idx if i < len(r) for a in r[i].split(",") if a.strip()}
  study = {l.strip() for p in pathlib.Path(".").glob("goatools_study_*.txt") for l in p.open(encoding="utf-8") if l.strip()}
  print(len(study), "study accessions,", len(study & outgroup), "from an outgroup column")
  PY
  ```

  Expected: `832 study accessions, 0 from an outgroup column`.

- `goatools_population.txt` still holds 483 lines, and
  `step3_final_consolidated_annotations.tsv` still holds 796 data rows — the
  rerun must not move those.

## Rollback

`git checkout -- step2/outputs/goatools step2/src/opxv_analysis/analyzer.py`.
