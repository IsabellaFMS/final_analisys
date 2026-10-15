# Step 4: host homology filtering

Removes from the viral proteome the proteins that share identical 12-residue
with the human proteome, so that later steps do not select epitopes
capable of cross-reacting with host sequences. The filtered FASTA is the input
of steps 5,6 and 7.

The k-mer comparison itself is not run here. It was executed separately against
human Uniprot entries, and this step consumes its output table: each row of
"sources/viral_hits.tsv" is one 12-mer shared between a viral protein and a
human one, with the position on both sides.

## Order of execution

Two things happen before "run_analysis.py":

1. **step 3**, which cleans the raw human and viral proteomes.
2. the external k-mer comparison, on the two cleaned FASTAs, producing
   "viral_hits.tsv".

"viral_hits.tsv" is a product of that comparison, so "run_analysis.py" cannot be
the entry point of the step on a fresh proteome.

## Inputs of run_analysis.py

In "sources/":

- "viral_hits.tsv", the hit table. Only the "query_seqname" column is read.
- "clean_viral_proteins.fasta", the proteome to be filtered, as cleaned by
  step 3. This is a local copy; the resolver prefers the original.

Two more inputs come from step 3 and drive the reincorporation described below.
They are **required**: the resolver looks for them under "../step3/outputs/"
and then under "sources/", and stops the step with an error if it finds
neither, rather than writing a proteome 45 records short without saying so.

- "removed_viral_proteins.csv", what step 3 removed and why.
- "removed_viral_proteins.fasta", the sequences it removed.

So step 3 has to have run. That was already true for "clean_viral_proteins.fasta"
in practice; now it is enforced.

## Usage

```bash
python run_analysis.py
```

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
line "of which identical to an excluded accession":

**33 of the 45 are identical to an accession this step excluded for host
homology.** Identical sequences share every 12-mer, so those 33 share the human
12-mers that got their twin excluded. They do not appear in "viral_hits.tsv"
only because the external comparison collapsed repeated sequences and reported
each hit once, through a single representative. Reincorporating them therefore
puts sequences with human homology back into "safe_sequences.fasta" — which the
191-record proteome produced before step 3 existed had always contained, and
which deduplication had hidden rather than caused. The remaining 12 pass the
filter on their own.

The exclusion list is unaffected: it is derived only from "viral_hits.tsv".
