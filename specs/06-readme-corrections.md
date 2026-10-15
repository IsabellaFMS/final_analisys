# 06 — README corrections

Covers TASK 4 of the brief. **Run last**, after specs 01–05: the READMEs must
describe the code as fixed, not as deposited.

Every item here is a factual error — the repository's own files contradict the
text — or an omission that leaves a manuscript claim unsupported.

## Divergences from the brief, verified first

The brief's line numbers were taken from a different revision of the root
README. Checked against the deposited file (91 lines):

| Brief says | Actually |
|---|---|
| line 29, background population | matches, line 29 |
| line 41, identical 12-residue stretches | matches, lines 40–41 |
| line 55, contiguous regions above the cutoff | matches, lines 54–56 |
| line 68, a data-flow table row | **no data-flow table exists** in this README |
| lines 112–122, tool versions "not recorded" | **no versions section exists**; the file ends at line 91 |

So two of the Task 4 items are additions, not edits. Treat them as such, and do
not invent the table the brief describes: add the minimum that removes the
error.

Also, as the brief instructs: **root line 3 ("Seven analysis steps…") stays as
is** — the manuscript's "six sequential modules" is what gets corrected. Note
that spec 04 adds step 8, so line 3 must say eight once that step exists. Make
that change as part of spec 04's documentation, and keep the brief's point: the
count follows the repository, and the manuscript follows the repository.

## `README.md` (root)

**Line 29** — replace

> The proteome of each virus is the background population.

with

> The background population is the set of 483 proteins carrying at least one GO
> term, pooled across the four proteomes.

The code uses `geneid2gos.keys()`, i.e. the contents of
`goatools_population.txt`, which has 483 lines.

**Lines 40–41** — replace

> …share identical 12-residue stretches with the human proteome…

with

> …contain a 12-residue window with 9 or more residues identical to a human
> 12-mer…

Only 1,078 of the 19,709 hit rows are exact matches; 13,995 are 9/12. Add the
parameters of the external k-mer run in the same place: window 12, minimum 9
identical residues, the human reference used and its release. The Maguire et al.
scripts are not redistributed, so this threshold is the only thing that makes
the filter reproducible — if a parameter cannot be recovered, say so explicitly
rather than omitting it (see [08-open-questions.md](08-open-questions.md) §2).

**Lines 54–56** — replace

> …extracts contiguous regions above the score cutoff, keeping those between 20
> and 30 residues.

with

> …masks runs of more than three consecutive sub-threshold residues, anchors
> each candidate at the highest remaining score and grows it symmetrically,
> keeping regions of 20 to 30 residues.

`step6/README.md` already describes this correctly. Align the root summary to
step 6, not the other way round — step 6 is the reference and spec 07 leaves it
untouched.

**New: the orthogroups file.** The brief wants the ambiguity in
`orthogroups.tsv` resolved. The deposited file has 339 rows
(`OG0000000`–`OG0000338`, 1,059 proteins); OrthoFinder's file of that name has
211. `step1/README.md` already explains the 128-id overlap; the root README says
nothing. Pick one:

- preferred: state in the root README that `step1/sources/orthogroups.tsv` is the
  concatenation of OrthoFinder's assigned and unassigned tables, pointing at
  `step1/README.md` for the overlap; or
- rename it to `orthogroups_merged.tsv` and update **every** reference —
  `step1/`, `step2/`'s R script, spec 02's enrichment input, spec 04's step 8,
  and all READMEs.

The rename is cleaner and riskier: the file name may be cited in the manuscript's
data-availability section. Default to the first option unless the author confirms
the name is not cited.

**New: tool versions.** Add a short "Tool versions" section listing what was
actually run, with the ones that are known:

- OrthoFinder v2.5.5, MCL inflation 1.7
- GOATOOLS v1.6.5
- eggNOG-mapper v2.1.12 (already in `step2/README.md`)
- InterProScan — **unresolved**, see §4 of the open questions; the deposited
  `interpro_orthopoxvirus.tsv` carries no version column, so neither
  `v5.75-106.019` nor `v5.75-106.0` can be confirmed from the repository
- DeepTMHMM, SignalP6, BepiPred-3.0, IEDB — record the version if it is
  recoverable; where it genuinely was not recorded, write "not recorded" rather
  than leaving it out

**New: three statements the brief requires, each verifiable:**

- 36 proteins (37–86 aa) are absent from
  `step2/outputs/tables/step3_final_consolidated_annotations.tsv` — neither
  InterProScan nor eggNOG-mapper returned anything for them. The table therefore
  has 796 rows, not 832, while the manuscript's denominator is 832, counting
  those 36 as unannotated.
- Pruning redundant parent GO terms upstream is statistically neutral given
  `propagate_counts=True`: intersecting the pruned sets and intersecting
  ancestor-closed sets both give 744 terms tested and FDR 0.0646.
- The run log reports 41 obsolete GO ids absent from the deposited `go-basic.obo`
  (release 2024-01-17); GOATOOLS skips them without error.

Confirm all three against the files before writing them. The 36-protein and
41-id figures must come out of the spec 02 rerun.

## `step2/README.md`

**Line 1** — `# Step 2: functional annotation and GO enrichment` →
`# Step 2: functional annotation and GO term distribution`. The reported result
is null, so "enrichment" in the title misleads.

**Line 11** — same replacement as root line 29, worded identically.

**Lines 57–58** — replace

> Of the 796 proteins, 241 carry at least one GO term: 150 with a biological
> process, 168 with a cellular component and 59 with a molecular function.

with

> Of the 796 proteins, 483 carry at least one GO term: 333 with a biological
> process, 230 with a cellular component and 323 with a molecular function.

These are stale numbers from an earlier run; the table in the same directory
gives the correct ones.

**Lines 69–70** — `Enrichment writes a result file only when at least one term
passes the FDR cutoff, so a missing file for a group means nothing was
significant for it.` This stops being true under spec 02, which writes every
tested term. Replace it with a description of the new behaviour.

**New: the two-level analysis.** Describe the protein level (exploratory; the
orthologues of one family are not independent observations) and the orthogroup
level (the manuscript's), with the union and intersection aggregation rules, the
850/744 terms tested, and the null result at FDR 0.05. List the
`outputs/enrichment_dual_level/` files, including the two partition lists.

## `step4/README.md`

- **Line 3**: the same `identical 12-residue` correction as root lines 40–41.
  Note that lines 51–65 of this README are about byte-identical duplicate
  *sequences*, which is a different claim and is correct — do not touch them.
- Record the external k-mer run's parameters, as above.
- **New**: name `outputs/reports/flagged_reintegrated_accessions.txt` (spec 03)
  where the text says the 33 are flagged, so the claim points at a file.

## `step6/README.md`

**No changes.** It is correct and is the reference for the root summary.

## `step8/README.md`

Written as part of spec 04. It must state the five rules in order, the parameter
values (0.2512, 0.2012, the truncation fraction, the accepted step 5 groups), the
reading chosen for each ambiguous rule, and the two selected antigens.

## Acceptance

Every number in every README traces to a file in the repository. Check at least:

```bash
grep -c . step2/outputs/goatools/goatools_population.txt            # 483
tail -n +2 step2/outputs/tables/step3_final_consolidated_annotations.tsv | grep -c .  # 796
grep -c . step4/sources/viral_hits.tsv                              # 19710 incl. header
grep -c . step4/outputs/reports/flagged_reintegrated_accessions.txt # 34 incl. header
awk -F'\t' 'NR>1{n++} END{print n}' step1/sources/orthogroups.tsv   # 339
```

- No README sentence describes behaviour the code no longer has.
- No README claim lacks a file that supports it — in particular the flagged
  accessions, the enrichment partitions and the step 8 selection trace.
- The root README's step count matches the number of steps present.
