# Step 3: proteome cleaning

Prepares the raw human and viral proteomes for the k-mer comparison and for the
filtering of step 4, by dropping the records that cannot take part in that
comparison: sequences that are duplicated, which would only multiply the same
k-mers and the same result rows, and sequences shorter than the k-mer size,
which produce no k-mer at all.

## Inputs

In `sources/`:

- `human.fasta`, the human proteome as downloaded from UniProt.
- `viral_proteins.fasta`, the RefSeq proteomes of the four orthopoxviruses.

## How cleaning works

Each proteome is treated independently, in this order:

1. **Too short.** A sequence shorter than `--kmer-size` is dropped, with reason
   `Short`. The default is 12, the k used in the comparison downstream. A
   sequence of exactly 12 residues yields one 12-mer and is kept.
2. **Duplicated.** Among the remaining records, sequences are grouped by exact
   match and only the **first occurrence in file order** is kept. The others are
   dropped with reason `Duplicate` and a `kept_instead` pointing at the survivor.


## Usage

```bash
python run_analysis.py
```

## Outputs

In `outputs/`:

- `fasta/clean_human_proteome.fasta` and `fasta/clean_viral_proteins.fasta`,
  the two cleaned proteomes.
- `reports/removed_human_proteins.csv` and `reports/removed_viral_proteins.csv`,
  one row per removed record with its reason, the accession kept in its place
  when the reason is `Duplicate`, and `source_index`, the position the record
  held in the raw file.
- `fasta/removed_viral_proteins.fasta`, the viral records that were removed.
- `reports/cleaning_summary.txt`, the run counts.

## Notes

- The k-mer comparison itself is not run here. This step produces the two FASTAs
  that it consumes; its output table is the input of step 4.

