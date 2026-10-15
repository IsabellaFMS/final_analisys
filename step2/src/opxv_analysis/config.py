from pathlib import Path

from .models import PipelinePaths


EGGNOG_COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]

GO_NAMESPACE_COLUMNS = {
    "biological_process": "GO_BP",
    "cellular_component": "GO_CC",
    "molecular_function": "GO_MF",
}

SUMMARY_GO_COLUMNS = ["Top 5 GOs (Combined)", "Top 10 GOs (Combined)"]
SUMMARY_REACTOME_COLUMNS = ["Top 5 Pathways (Reactome)", "Top 10 Pathways (Reactome)"]

SPECIES_MAP = {
    "cowpox_refseq_Brighton_Red_proteins": "CPXV",
    "mpox_refseq_Rivers_State_proteins": "MPXV",
    "vaccinia_WR_refseq_proteins": "VACV",
    "varv_refseq_proteins": "VARV",
}

# Columns of step1/sources/orthogroups.tsv that are not part of the analysis.
OUTGROUP_COLUMNS = ("Humanalphaherpesvirus1_proteins", "swinepox_proteins")

# Abbreviations accepted inside the txt conservation file names.
CONSERVATION_SPECIES = tuple(SPECIES_MAP.values())


def default_paths(root: Path | None = None) -> PipelinePaths:
    root = root or Path(__file__).resolve().parents[2]
    outputs = root / "outputs"
    return PipelinePaths(
        root=root,
        sources=root / "sources",
        outputs=outputs,
        tables=outputs / "tables",
        reports=outputs / "reports",
        texts=outputs / "texts",
        visualizations=outputs / "visualizations",
        goatools=outputs / "goatools",
        enrichment=outputs / "enrichment_dual_level",
    )
