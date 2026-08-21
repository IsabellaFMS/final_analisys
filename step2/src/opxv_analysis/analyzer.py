import re
from collections import defaultdict

import pandas as pd

from .config import GO_NAMESPACE_COLUMNS
from .loader import sanitize_group_name
from .models import ComparisonStats, PipelinePaths


def locate_column(columns: list[str], *terms: str) -> str:
    for column in columns:
        column_norm = column.lower()
        if all(term in column_norm for term in terms):
            return column
    raise KeyError(f"Column not found for search terms: {terms}")


def combine_go_terms(interpro_value: str, eggnog_value: str) -> str | None:
    values = set()
    for part in str(interpro_value).split("|"):
        part = part.strip()
        if part and part != "-" and part.startswith("GO:"):
            values.add(part)
    for part in str(eggnog_value).split(","):
        part = part.strip()
        if part and part != "-" and part.startswith("GO:"):
            values.add(part)
    return "|".join(sorted(values)) if values else None


def merge_annotations(eggnog_df: pd.DataFrame, interpro_df: pd.DataFrame) -> pd.DataFrame:
    
    interpro = interpro_df.copy()
    selected = interpro[[
        "Protein accession",
        locate_column(list(interpro.columns), "interpro annotations", "accession"),
        locate_column(list(interpro.columns), "interpro annotations", "description"),
        locate_column(list(interpro.columns), "go annotations"),
        locate_column(list(interpro.columns), "pathways annotations"),
    ]].copy()
    selected.columns = [
        "Protein_Accession", "Interpro_Accession", "Interpro_Description",
        "GO_Annotations", "Pathways_Annotations",
    ]

    eggnog = eggnog_df.rename(columns={"query": "Protein_Accession"})[["Protein_Accession", "Description", "GOs"]]
    eggnog = eggnog.drop_duplicates(subset=["Protein_Accession"]).copy()
    merged = pd.merge(selected, eggnog, on="Protein_Accession", how="outer")

    
    merged["Combined_GO_Annotations"] = merged.apply(lambda row: combine_go_terms(row.get("GO_Annotations"), row.get("GOs")), axis=1)
    final = merged[[
        "Protein_Accession", "Interpro_Accession", "Interpro_Description", "Description",
        "Combined_GO_Annotations", "Pathways_Annotations",
    ]].copy()
    final.rename(columns={"Description": "eggNOG_Description"}, inplace=True)
    return final


def strip_go_source_tag(term: str) -> str:
    # InterPro reports its terms as GO:0004222(InterPro), so the source tag has
    # to come off before the identifier can be found in the ontology.
    match = re.match(r"^(GO:\d+)", term.strip())
    return match.group(1) if match else term.strip()


def refine_go_list(go_string: str, godag) -> str | None:
    
    if pd.isna(go_string) or not str(go_string).strip():
        return None
    original = {
        go_id
        for go_id in (strip_go_source_tag(term) for term in str(go_string).split("|"))
        if go_id in godag
    }
    if not original:
        return None
    to_remove = set()
    for go_id in original:
        to_remove.update(original.intersection(godag[go_id].get_all_parents()))
    refined = original - to_remove
    return "|".join(sorted(refined)) if refined else None


def annotate_reactome(pathway_string: str, reactome_map: dict[str, str]) -> str | None:
    if pd.isna(pathway_string):
        return None
    descriptions = []
    for part in re.split(r"[|,]", str(pathway_string)):
        part = part.strip()
        if ":" not in part:
            continue
        db, path_id = part.split(":", 1)
        if db.lower() == "reactome" and path_id in reactome_map:
            descriptions.append(reactome_map[path_id])
    return " ; ".join(sorted(set(descriptions))) if descriptions else None


def refine_annotations(df: pd.DataFrame, godag, reactome_map: dict[str, str]) -> pd.DataFrame:
    refined = df.copy()
    refined["Refined_GO_Annotations"] = refined["Combined_GO_Annotations"].apply(lambda value: refine_go_list(value, godag))
    refined["Pathways_Description"] = refined["Pathways_Annotations"].apply(lambda value: annotate_reactome(value, reactome_map))
    return refined


def extract_unique_values(series: pd.Series) -> str:
    values = set()
    for item in series.dropna():
        if item and item != "-":
            for part in re.split(r"\s*;\s*|\s*\|\s*", str(item)):
                if part:
                    values.add(part.strip())
    return " ; ".join(sorted(values)) if values else "-"


def split_go_by_namespace(go_string: str, godag) -> dict[str, str]:
    # One column per GO namespace, each term carrying its own name.
    buckets: dict[str, list[str]] = {column: [] for column in GO_NAMESPACE_COLUMNS.values()}
    if go_string and go_string != "-":
        for go_id in str(go_string).split(" ; "):
            term = godag.get(go_id.strip())
            if term is None:
                continue
            column = GO_NAMESPACE_COLUMNS.get(term.namespace)
            if column:
                buckets[column].append(f"{go_id.strip()} ({term.name})")
    return {column: " ; ".join(values) if values else "-" for column, values in buckets.items()}


def count_terms_outside_dag(series: pd.Series, godag) -> int:
    # Terms the OBO does not know have no namespace, so they fit in no column.
    missing = 0
    for value in series:
        if not value or value == "-":
            continue
        missing += sum(1 for go_id in str(value).split(" ; ") if go_id.strip() not in godag)
    return missing


def consolidate_annotations(df: pd.DataFrame, godag) -> pd.DataFrame:
    grouped = df.fillna("").groupby("Protein_Accession").agg({
        "Interpro_Accession": extract_unique_values,
        "Interpro_Description": extract_unique_values,
        "eggNOG_Description": extract_unique_values,
        "Refined_GO_Annotations": extract_unique_values,
    }).reset_index()

    dropped = count_terms_outside_dag(grouped["Refined_GO_Annotations"], godag)
    if dropped:
        print(f"Warning: {dropped} GO term(s) absent from the OBO were dropped from the namespace columns.", flush=True)

    namespaces = grouped["Refined_GO_Annotations"].apply(lambda value: pd.Series(split_go_by_namespace(value, godag)))
    grouped = pd.concat([grouped, namespaces], axis=1)
    return grouped[["Protein_Accession", "eggNOG_Description", "Interpro_Accession", "Interpro_Description", *GO_NAMESPACE_COLUMNS.values()]]


def compare_annotation_coverage(eggnog_df: pd.DataFrame, interpro_df: pd.DataFrame) -> ComparisonStats:
    # Count how many proteins each tool annotates, and how many both cover.
    
    interpro = interpro_df.copy()
    interpro["Score"] = pd.to_numeric(interpro["Score"], errors="coerce")


    no_filter = interpro["Analysis"].isin(["ProSitePatterns", "Coils", "MobiDB Lite", "SMART"])
    evalue_filter = interpro["Analysis"].isin(["Pfam", "TIGRFAMs", "SUPERFAMILY", "PANTHER", "ProSiteProfiles"]) & (interpro["Score"] < 1e-5)
    bitscore_filter = interpro["Analysis"].isin(["Gene3D", "CDD"]) & (interpro["Score"] > 50)
    interpro_ids = set(interpro[no_filter | evalue_filter | bitscore_filter]["Protein accession"].dropna().str.strip().unique())
    eggnog_ids = set(eggnog_df["query"].dropna().str.strip().unique())
    union_total = len(eggnog_ids.union(interpro_ids))
    both = len(eggnog_ids.intersection(interpro_ids))
    return ComparisonStats(len(eggnog_ids), len(interpro_ids), union_total, both, len(interpro_ids - eggnog_ids), len(eggnog_ids - interpro_ids))


def build_goatools_inputs(eggnog_df: pd.DataFrame, interpro_df: pd.DataFrame, conservation_df: pd.DataFrame, paths: PipelinePaths) -> list[str]:

    eggnog = eggnog_df.rename(columns={"query": "protein_id"}).copy()
    interpro = interpro_df.rename(columns={"Protein accession": "protein_id"}).copy()
    conservation = conservation_df.rename(columns={"Protein_ID": "protein_id", "Conservation_Group": "conservation_group"}).copy()

    eggnog_go = eggnog[eggnog["GOs"].fillna("-") != "-"].copy()
    interpro_go = interpro[interpro["GO annotations"].fillna("-") != "-"].copy()

    
    go_map: dict[str, set[str]] = defaultdict(set)
    for _, row in eggnog_go.iterrows():
        go_map[row["protein_id"]].update([part.strip() for part in str(row["GOs"]).split(",") if part.strip()])
    for _, row in interpro_go.iterrows():
        valid = [term.split("(")[0] for term in str(row["GO annotations"]).split("|") if term.startswith("GO:")]
        if valid:
            go_map[row["protein_id"]].update(valid)

    population_file = paths.goatools / "goatools_population.txt"
    population_file.write_text("".join(f"{protein}\t{';'.join(sorted(values))}\n" for protein, values in go_map.items()), encoding="utf-8")

    created = [population_file.name]
    unique_groups = conservation.drop_duplicates(subset=["protein_id"])["conservation_group"].unique()
    for group in unique_groups:
        group_file = paths.goatools / f"goatools_study_{sanitize_group_name(group)}.txt"
        proteins = conservation[conservation["conservation_group"] == group]["protein_id"].tolist()
        group_file.write_text("\n".join(proteins), encoding="utf-8")
        created.append(group_file.name)
    return created


def format_top_items(series: pd.Series, top_n: int = 5) -> str:
    cleaned = series.dropna()
    cleaned = cleaned[~cleaned.isin(["-", "S", ""])].astype(str).str.strip()
    if cleaned.empty:
        return "None found"
    return ", ".join(f"{name} ({count})" for name, count in cleaned.value_counts().head(top_n).items())


def build_group_summary(eggnog_df: pd.DataFrame, interpro_df: pd.DataFrame, conservation_df: pd.DataFrame, total_proteins: int | None = None) -> tuple[pd.DataFrame, dict[str, float | int]]:
    eggnog = eggnog_df.rename(columns={"query": "protein_id"}).copy()
    interpro = interpro_df.rename(columns={"Protein accession": "protein_id"}).copy()
    conservation = conservation_df.rename(columns={"Protein_ID": "protein_id", "Conservation_Group": "conservation_group"}).copy()
    conservation_unique = conservation.drop_duplicates(subset=["protein_id"]).copy()

    eggnog_go = eggnog[eggnog["GOs"].fillna("-") != "-"].copy()
    interpro_go = interpro[interpro["GO annotations"].fillna("-") != "-"].copy()
    proteins_eggnog_go = set(eggnog_go["protein_id"].unique())
    proteins_interpro_go = set(interpro_go["protein_id"].unique())
    all_with_go = proteins_eggnog_go.union(proteins_interpro_go)
    total_proteins = total_proteins or conservation_unique["protein_id"].nunique()

    rows = []
    for group in sorted(conservation_unique["conservation_group"].unique()):
        protein_ids = set(conservation_unique[conservation_unique["conservation_group"] == group]["protein_id"])
        interpro_group = interpro[interpro["protein_id"].isin(protein_ids)]
        eggnog_group = eggnog[eggnog["protein_id"].isin(protein_ids)]
        interpro_gos = interpro_group["GO annotations"].str.split("|").explode().str.split("(").str[0]
        eggnog_gos = eggnog_group["GOs"].str.split(",").explode()
        annotated = protein_ids.intersection(all_with_go)
        rows.append({
            "Group": group,
            "Total Proteins": len(protein_ids),
            "With GO (%)": f"{len(annotated)} ({len(annotated)/len(protein_ids):.1%})",
            "Top 5 COG Categories": format_top_items(eggnog_group["COG_category"], 5),
            "Top 5 Domains (InterPro)": format_top_items(interpro_group["Signature description"], 5),
            "Top 10 Domains (InterPro)": format_top_items(interpro_group["Signature description"], 10),
            "Top 5 GOs (Combined)": format_top_items(pd.concat([interpro_gos, eggnog_gos]), 5),
            "Top 10 GOs (Combined)": format_top_items(pd.concat([interpro_gos, eggnog_gos]), 10),
        })

    stats = {
        "total_proteins": total_proteins,
        "annotated_proteins": len(all_with_go),
        "unannotated_proteins": total_proteins - len(all_with_go),
        "eggnog_go": len(proteins_eggnog_go),
        "interpro_go": len(proteins_interpro_go),
        "both_go": len(proteins_eggnog_go.intersection(proteins_interpro_go)),
        "eggnog_only_go": len(proteins_eggnog_go - proteins_interpro_go),
        "interpro_only_go": len(proteins_interpro_go - proteins_eggnog_go),
    }
    return pd.DataFrame(rows), stats
