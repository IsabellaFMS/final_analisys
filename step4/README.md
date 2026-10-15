# Step 4: host homology filtering

Removes from the viral proteome the proteins that contain a 12-residue window
with 9 or more residues identical to a human 12-mer, so that later steps do not
select epitopes capable of cross-reacting with host sequences. The filtered
FASTA is the input of steps 5, 6 and 7.

The k-mer comparison itself is not run here. It was executed separately against
human UniProt entries, with a window of 12 residues and an acceptance threshold
of 9 identical residues out of 12, and this step consumes its output table: each
row of `sources/viral_hits.tsv` is one 12-residue window shared between a viral
protein and a human one, with the position on both sides.

Of the 19,709 rows, 1,078 are exact 12/12 matches and 13,995 match at 9/12,
over 641 distinct viral proteins. The threshold matters: the comparison is not
restricted to identical stretches, and describing it that way would understate
how much of the proteome the filter removes. The version of the external
comparison scripts and the release of the human reference were not recorded,
so the deposited hit table and this threshold are what makes the filter
reproducible.

## Order of execution

Two things happen before `run_analysis.py`:

1. **step 3**, which cleans the raw human and viral proteomes.
2. the external k-mer comparison, on the two cleaned FASTAs, producing
   `viral_hits.tsv`.

`viral_hits.tsv` is a product of that comparison, so `run_analysis.py` cannot be
the entry point of the step on a fresh proteome.

## Inputs of run_analysis.py

In `sources/`:

- `viral_hits.tsv`, the hit table. Only the `query_seqname` column is read.
- `clean_viral_proteins.fasta`, the proteome to be filtered, as cleaned by
  step 3. This is a local copy; the resolver prefers the original.

Two more inputs come from step 3 and drive the reincorporation described below.
They are **required**: the resolver looks for them under `../step3/outputs/`
and then under `sources/`, and stops the step with an error if it finds
neither, rather than writing a proteome 45 records short without saying so.

- `removed_viral_proteins.csv`, what step 3 removed and why.
- `removed_viral_proteins.fasta`, the sequences it removed.

So step 3 has to have run. That was already true for `clean_viral_proteins.fasta`
in practice; now it is enforced.

## Usage

```bash
python run_analysis.py
```

## Outputs

- `outputs/reports/exclusion_ids.txt`, the 641 accessions with at least one hit
- `outputs/fasta/safe_sequences.fasta`, the filtered proteome, 191 records
- `outputs/reports/flagged_reintegrated_accessions.txt`, the 33 reincorporated
  duplicates identical to an excluded accession, as a tab-separated table with
  `protein_id`, `kept_instead`, `pathogen` and `source_index`

## Reincorporation of the step 3 duplicates

Step 3 removes sequences that are byte-identical to one it keeps, which takes
the filtered proteome from 191 to 146 records. Those 45 accessions carry their
own species and their own conservation class, and the sequence content they
represent is not lost by deduplication precisely because it is identical, so
this step puts all 45 back and writes 191 records again. Each is reinserted at
the position it held in the raw proteome, which keeps the output identical, in
content and in order, to the FASTA the external predictors of steps 5 to 7 were
run on.

One caveat travels with the 45, and the execution summary counts it on the
line `of which identical to an excluded accession`:

**33 of the 45 are identical to an accession this step excluded for host
homology**, and they are listed in
`outputs/reports/flagged_reintegrated_accessions.txt`, one row each with the
accession kept in their place, the pathogen and the position they held in the
raw proteome. Identical sequences share every 12-mer, so those 33 share the human
12-mers that got their twin excluded. They do not appear in `viral_hits.tsv`
only because the external comparison collapsed repeated sequences and reported
each hit once, through a single representative. Reincorporating them therefore
puts sequences with human homology back into `safe_sequences.fasta` — which the
191-record proteome produced before step 3 existed had always contained, and
which deduplication had hidden rather than caused. The remaining 12 pass the
filter on their own.

The exclusion list is unaffected: it is derived only from `viral_hits.tsv`.
