# Step 8: final antigen selection

Applies the final selection rules to the step 5, 6 and 7 tables and reports the
antigens that survive all five. This is the step that produces the shortlist; it
runs no predictor and introduces no new data.

## Inputs

Each one is resolved from its own step, so nothing has to be copied into
`sources/`:

- `../step5/outputs/tables/comparative_results_signalp6_deep_final.csv`, the
  SignalP6 / DeepTMHMM consensus and the conservation label
- `../step6/run_analysis.py`, used to extract the epitopes at the two
  stringencies this step needs (see below)
- `../step7/outputs/mhci_final.csv` and `../step7/outputs/mhcii_final.csv`, for
  allelic coverage
- `../step1/sources/orthogroups.tsv`, for orthogroup membership and the species
  present in each orthogroup. The `Humanalphaherpesvirus1_proteins` and
  `swinepox_proteins` columns are excluded, as in step 1
- `../step2/sources/all_refseq_proteins.fasta`, for sequence lengths and
  descriptions

Lengths come from the full 832-protein proteome rather than from the filtered
one, because the orthogroup consensus length has to be computed over every
member. Taking it over the filtered proteome would let a truncated orthologue
look full-length against a depleted orthogroup.

## The two B-cell stringencies

The deposited `step6/outputs/tables/epitopes.csv` is the 0.1512 run. This step
needs 0.2512 (the selection criterion) and 0.2012 (the fallback reporting
stringency of rule 5), so it calls step 6 with its own `--cutoff` and
`--output-dir`:

```bash
python ../step6/run_analysis.py --cutoff 0.2512 --output-dir outputs/bcell/0.2512
```

The extraction logic is step 6's, unmodified, and the deposited 0.1512 table is
not touched. The tables land under `outputs/bcell/<cutoff>/tables/epitopes.csv`,
and the per-protein FASTA files step 6 writes alongside them land in the same
tree — that is step 6's output contract, not this step's.

Existing tables are reused; `--refresh` re-runs the extraction. The three runs
reproduce 473 epitopes across 168 proteins at 0.1512, 168 across 97 at 0.2012
and 40 across 28 at 0.2512.

## The rules

Applied in order, to the 191 proteins of the step 4 proteome:

1. **B-cell evidence.** Keep proteins with at least one epitope at cutoff
   0.2512. 28 proteins.
2. **Conservation.** Keep proteins whose orthogroup holds all four OPXV species.
   16 proteins.
3. **Surface accessibility.** Keep `Group 1: High Confidence Consensus
   (SP+/THMM+)`, the agreement of SignalP6 and DeepTMHMM. 3 proteins. Groups 2
   and 3 are single-predictor positives and are not treated as a consensus; the
   accepted set is the `ACCESSIBLE_GROUPS` constant in `src/antigen_selection/config.py`.
4. **One representative per orthogroup.** Where more than one orthologue passes,
   keep the richest epitope profile: most epitopes at 0.2512, ties broken by the
   combined count of distinct MHC-I and MHC-II alleles, then by accession so the
   result does not depend on input order. 2 orthogroups.
5. **Truncation.** A representative shorter than half its orthogroup consensus
   length, or described as partial, does not represent its orthogroup. Where it
   is the only member meeting the 0.2512 criterion, the orthogroup is
   represented instead by its full-length orthologue — surface accessible, with
   at least one epitope at 0.2012 — ranked by the same criterion as rule 4, and
   its epitopes are reported at 0.2012.

The truncation threshold is `TRUNCATION_RATIO = 0.5`. In the one orthogroup
where the rule fires, the truncated member sits at 0.35 of the consensus and the
next shortest surviving member at 1.00, so the outcome is not sensitive to the
exact value within 0.4 to 0.5.

## Usage

```bash
python run_analysis.py
```

## Outputs

Under `outputs/`:

- `tables/final_antigens.tsv`, the selected antigens with their orthogroup,
  lengths, the stringency their epitopes are reported at, epitope count, allelic
  coverage, accessibility group, and the accession each one replaces, if any
- `tables/selection_trace.tsv`, one row per protein that entered rule 1 plus any
  substitute, with the verdict that decided it and the reason. This is the audit
  trail for the selection
- `tables/selected_epitopes.tsv`, the epitopes of the selected antigens with the
  step 6 columns unchanged, plus the cutoff they came from
- `reports/selection_summary.txt`, the counts surviving each rule, the
  parameters used and the orthogroups where rules 4 and 5 fired
- `bcell/<cutoff>/`, the step 6 extraction at each stringency

## Result

Two antigens:

| Protein | Species | Orthogroup | Epitopes | Stringency |
|---|---|---|---|---|
| `NP_619793.1`, CrmB (CPXV005) | CPXV | `OG0000003` | 3 | 0.2012 |
| `YP_010377182.1`, chemokine-binding protein | MPXV | `OG0000008` | 4 | 0.2512 |

Rule 5 fires for `OG0000003`: its only member with an epitope at 0.2512 is
`YP_233097.1`, 122 aa against a consensus of 349 aa, so the orthogroup is
represented by the full-length `NP_619793.1` and its epitopes are reported at
0.2012. `NP_619793.1` is chosen over the other full-length accessible member,
`YP_010377181.1`, on epitope count at 0.2012 (3 against 2).

Both selected accessions are reincorporated duplicates listed in
`../step4/outputs/reports/flagged_reintegrated_accessions.txt`: each is identical
in sequence to an accession the host-homology filter excluded, and therefore
carries the same human 12-mer matches.
