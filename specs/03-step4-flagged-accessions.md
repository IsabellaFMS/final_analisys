# 03 — Step 4: persist the flagged accessions

Covers TASK 3 of the brief. Independent of specs 01, 02 and 04; can run in
parallel with them.

## Problem

`step4/README.md` states that 33 of the 45 reincorporated duplicates are
identical to accessions excluded for host homology, and that they are flagged.
The manuscript's Results section says the same, and says they are flagged *in the
deposited data*.

The computation exists and is correct:

- `step4/src/viral_filtering/analyzer.py:88` builds
  `DuplicateVerdict(..., twin_excluded=kept_instead in exclusion_ids, ...)`;
- `ReincorporationSummary.twin_excluded_count`
  (`step4/src/viral_filtering/models.py:80`) counts them;
- `reporting.py:31` prints `of which identical to an excluded accession: 33`.

But the list is never written. The verdicts live only in the returned
`StepResult`, which `run_analysis.py` discards after printing the summary. So
nothing in the deposit is in fact flagged, and the manuscript's sentence has no
support a reviewer can check.

## Do this

Write `step4/outputs/reports/flagged_reintegrated_accessions.txt` from the
existing `ReincorporationSummary.verdicts` — no recomputation, no new logic in
`analyzer.py`.

- **`writer.py`**: add `write_flagged_verdicts(verdicts, output_path)` next to
  `write_exclusion_ids`, following the same shape (`_prepare_parent`, explicit
  `encoding="utf-8"`). Tab-separated with a header row:

  ```
  protein_id	kept_instead	pathogen	source_index
  ```

  One row per verdict with `twin_excluded` true, in ascending `source_index`
  order so the file is stable across runs and diffable. Write the file even when
  no verdict qualifies — a header-only file is a result; a missing file is
  ambiguous.

  The extension is `.txt` because the brief names that path, and
  `exclusion_ids.txt` in the same directory is also a text table. Keep the
  header so the columns are self-describing.

- **`services.py`**: call it from `_filter_and_write`, immediately after
  `write_fasta_records`, so both the `run` and the `filter` workflows produce it.
  The output path belongs in `CommandPaths` alongside `exclusion_output`, with
  default `outputs/reports/flagged_reintegrated_accessions.txt`, resolved in
  `config.resolve_command_paths` the way the other outputs are, and exposed as
  `--flagged-output` on the `run` and `filter` subparsers in `run_analysis.py`
  for symmetry with its siblings.

- **`reporting.py`**: no change to the counts, but the summary should name the
  file it just wrote, so the operator knows where the 33 went. One line, in the
  existing style.

## Acceptance

```bash
cd step4 && python run_analysis.py run
grep -c . outputs/reports/flagged_reintegrated_accessions.txt   # 34 (header + 33)
grep -E "YP_010377182\.1|NP_619793\.1" outputs/reports/flagged_reintegrated_accessions.txt
```

- 33 data rows.
- `YP_010377182.1` and `NP_619793.1` are both present. These are the two antigens
  the manuscript prioritises and the two that spec 04 must select, so a reviewer
  will look for exactly these.
- The printed summary still reports `Duplicates reincorporated: 45` and
  `of which identical to an excluded accession: 33`, and the row count of the new
  file equals the second number.
- `outputs/fasta/safe_sequences.fasta` and `outputs/reports/exclusion_ids.txt`
  are byte-identical to the baseline: this spec adds an output, it does not
  change the filtering.

## Rollback

`git checkout -- step4/` and delete the new file.
