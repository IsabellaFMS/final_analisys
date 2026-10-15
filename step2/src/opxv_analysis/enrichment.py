import argparse

from goatools.associations import read_associations
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag


def run_go_enrichment(study_file: str, population_file: str, obo_file: str, output_file: str) -> None:

    
    print(f"--- Starting analysis for: {study_file} ---")
    obodag = GODag(obo_file)

    study_gene_ids = [line.strip() for line in open(study_file, encoding="utf-8", errors="replace") if line.strip()]
    geneid2gos = read_associations(population_file, "id2gos")
    if not study_gene_ids:
        print("WARNING: the study file is empty. Skipping analysis.")
        return

    # Fisher's exact test with Benjamini-Hochberg correction at FDR 0.05, the configuration reported in the methods. 
    # propagate_counts pulls the counts of each term up to its ancestors, so a parent term is not called depleted only because the annotation sits on a more specific child.
    goeaobj = GOEnrichmentStudy(
        geneid2gos.keys(),
        geneid2gos,
        obodag,
        propagate_counts=True,
        alpha=0.05,
        methods=["fdr_bh"],
    )
    results = goeaobj.run_study(study_gene_ids)
    significant = [result for result in results if result.p_fdr_bh < 0.05]
    if significant:
        goeaobj.wr_tsv(output_file, significant)
        print(f"Results saved to: {output_file}")
    else:
        print("Analysis completed. No significantly enriched GO term was found.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run GO enrichment analysis with GOATOOLS.")
    parser.add_argument("study_file")
    parser.add_argument("population_file")
    parser.add_argument("obo_file")
    parser.add_argument("output_file")
    return parser
