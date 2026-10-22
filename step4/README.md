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
  step 3.

## Usage

```bash
python run_analysis.py
```

## Outputs

- `outputs/reports/exclusion_ids.txt`
- `outputs/fasta/safe_sequences.fasta`
- `outputs/reports/flagged_reintegrated_accessions.txt`

