# 08 — Open questions for the author

Five points the specs cannot settle from the repository alone. Each one affects
text that would go into a README or the manuscript, so each needs an answer
before spec 06 is finalised. None of them should be resolved by guessing.

## 1. `YP_233097.1` — species label

The brief (TASK 6, acceptance) describes the truncated CrmB orthologue as
"`YP_233097.1` (VARV, 122 aa against a 349–355 aa consensus)".

Verified against the deposit:

- the length is right: 122 aa, against members of 349, 349, 349, 355, 355, 320,
  202 and 122 in `OG0000003`;
- the species is **not**: `YP_233097.1` sits in the
  `vaccinia_WR_refseq_proteins` column of `step1/sources/orthogroups.tsv`, and
  `step6/outputs/fasta/` names it
  `YP_233097.1_TNF-alpha-receptor-like_protein_Vaccinia_virus.fasta`. The VARV
  member of that orthogroup is `NP_042240.1`, at 349 aa.

**Question.** Is the VARV label an error in the brief only, or does it also
appear in the manuscript's Final Antigen Selection section? If it is in the
manuscript, it needs correcting there — the selection logic is unaffected.

## 2. The external k-mer run

Specs 06 (root lines 40–41, `step4/README.md`) require the parameters of the
host-homology comparison, because the Maguire et al. scripts are not
redistributed and the deposit therefore cannot be reproduced without them.

Recoverable from `step4/sources/viral_hits.tsv`: 19,709 rows, 1,078 at 12/12
identical residues, 13,995 at 9/12, 641 distinct viral proteins with at least
one hit — so the window is 12 and the acceptance threshold is 9.

**Not** recoverable from the deposit:

- which human reference was used and its release (`step3/sources/human.fasta` is
  present, but its provenance and date are not recorded anywhere);
- the version of the Maguire et al. scripts, and any non-default setting;
- when the run was performed.

**Question.** Supply these, or state in both READMEs that they were not
recorded. The second is acceptable and honest; silence is not.

## 3. Rule 3 of step 8 — what "the consensus" admits

The manuscript says surface accessibility "per the SignalP6 / DeepTMHMM
consensus". Step 5 produces four groups:

| Group | Meaning | Proteins |
|---|---|---|
| 1 | High confidence consensus (SP+/THMM+) | 25 |
| 2 | SignalP only (SP+/THMM−) | 3 |
| 3 | DeepTMHMM only (SP−/THMM+) | 7 |
| 4 | Negative consensus (SP−/THMM−) | 156 |

(191 proteins in total, the step 4 output.)

Spec 04 implements Group 1 only, which is the reading "consensus" most naturally
supports and which both target antigens satisfy. But Groups 2 and 3 are
single-predictor positives, and a reader could take "consensus" to mean the
consolidated call rather than agreement between the two predictors.

**Question.** Confirm Group 1 only. If Groups 2 and 3 are also admitted, spec 04
must be re-run and the candidate set after rule 3 will be larger — the final two
may still survive rules 4 and 5, but that has to be demonstrated, not assumed.

## 4. InterProScan version

The brief says the README reads `v5.75-106.019` and the manuscript
`v5.75-106.0`, and asks which is right.

Neither can be confirmed from the repository: the deposited root README has no
versions section at all (it ends at line 91), and
`step2/sources/interpro_orthopoxvirus.tsv` carries no version column — its
columns are accession, MD5, length, analysis, signature, description, start,
stop, score, status, date, InterPro accession and description, GO and pathway
annotations. The analysis date in the file is `21-04-2025`.

**Question.** Supply the version from the InterProScan run log or the job
record. `v5.75-106.0` is a plausible release string and `v5.75-106.019` is not,
but that is an inference, not a fact from the deposit.

## 5. Rule 4 of step 8 — "more conserved B-cell epitopes"

The tie-break reads "more conserved B-cell epitopes, ties broken by combined MHC
class I and class II allelic coverage". "Conserved" is ambiguous:

- (a) epitopes belonging to a protein whose conservation class is `Conserved` —
  degenerate, because rule 2 has already restricted the set to conserved
  orthogroups, so every surviving candidate qualifies and the criterion decides
  nothing;
- (b) epitopes whose **sequence** is conserved across the orthologues of the
  orthogroup — meaningful, but it requires aligning the orthologues, which no
  step in the deposit does.

Spec 04 implements the count of epitopes at the 0.2512 cutoff, with allele
coverage as the tie-break, and records the choice. That is reading (a) reduced to
its only non-degenerate content.

**Question.** If the manuscript means (b), step 8 needs an alignment step and the
tie-break must be recomputed. Say which, and if (b), supply the alignment tool
and parameters used for the published result.
