from __future__ import annotations

from typing import Iterable

import pandas as pd


def _unique_join(values: Iterable[object]) -> str:
    clean = sorted({str(value).strip() for value in values if pd.notna(value) and str(value).strip()})
    return "; ".join(clean)


def _numeric_join(values: Iterable[object]) -> str:
    numbers = sorted({float(value) for value in values if pd.notna(value)})
    rendered = [str(int(value)) if value.is_integer() else f"{value:g}" for value in numbers]
    return "; ".join(rendered)


def concatenate_tables(frames: list[pd.DataFrame], *, deduplicate: bool = False) -> pd.DataFrame:
    combined = pd.concat(frames, ignore_index=True)
    if deduplicate:
        combined = combined.drop_duplicates().reset_index(drop=True)
    return combined


def filter_hits_by_rank(frame: pd.DataFrame, cutoff: float, rank_column: str = "rank") -> pd.DataFrame:
    ranked = frame.copy()
    ranked[rank_column] = pd.to_numeric(ranked[rank_column], errors="coerce")
    ranked = ranked[ranked[rank_column] <= cutoff].copy()
    return ranked.sort_values(by=rank_column, ascending=True).reset_index(drop=True)


def add_protein_ids(frame: pd.DataFrame, id_map: dict[int, str]) -> pd.DataFrame:
    enriched = frame.copy()
    enriched["protein_id"] = enriched["seq_num"].map(id_map)
    columns = ["protein_id"] + [column for column in enriched.columns if column != "protein_id"]
    return enriched[columns]


def analyze_mhci_promiscuity(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for peptide, group in frame.groupby("peptide"):
        alleles = sorted(group["allele"].dropna().astype(str).unique())
        rows.append(
            {
                "peptide": peptide,
                "allele_count": len(alleles),
                "best_rank": pd.to_numeric(group["rank"], errors="coerce").min(),
                "bound_alleles": "; ".join(alleles),
                "source_proteins": _unique_join(group["protein_id"]),
                "start": _numeric_join(group["start"]),
                "end": _numeric_join(group["end"]),
            }
        )
    result = pd.DataFrame(rows)
    return result.sort_values(by=["allele_count", "best_rank"], ascending=[False, True]).reset_index(drop=True)


def analyze_mhcii_promiscuity(frame: pd.DataFrame) -> pd.DataFrame:
    analyzed = frame.copy()
    analyzed["rank"] = pd.to_numeric(analyzed["rank"], errors="coerce")
    grouping_columns = ["protein_id", "core_peptide"]
    counts = analyzed.groupby(grouping_columns)["allele"].nunique()
    analyzed["core_allele_count"] = analyzed.set_index(grouping_columns).index.map(counts)
    return analyzed.sort_values(by=["core_allele_count", "rank"], ascending=[False, True]).reset_index(drop=True)


def build_mhci_final_output(frame: pd.DataFrame, classification_map: dict[str, str]) -> pd.DataFrame:
    output = frame.copy()
    output["score"] = pd.to_numeric(output["score"], errors="coerce")
    output["rank"] = pd.to_numeric(output["rank"], errors="coerce")

    grouped_rows: list[dict[str, object]] = []
    for (protein_id, allele, peptide), group in output.groupby(
        ["protein_id", "allele", "peptide"],
        dropna=False,
    ):
        grouped_rows.append(
            {
                "Protein_ID": protein_id,
                "Epitope": peptide,
                "Allele": allele,
                "Elution Score": group["score"].max(),
                "Rank": group["rank"].min(),
                "start": group["start"].dropna().iloc[0] if group["start"].notna().any() else pd.NA,
                "end": group["end"].dropna().iloc[0] if group["end"].notna().any() else pd.NA,
            }
        )

    output = pd.DataFrame(grouped_rows)
    if output.empty:
        columns = [
            "Protein_ID",
            "Allele",
            "Epitope",
            "Allele Count",
            "Elution Score",
            "Rank",
            "Alleles",
            "Classification (conservation)",
            "start",
            "end",
        ]
        return pd.DataFrame(columns=columns)

    grouping_columns = ["Protein_ID", "Epitope"]
    allele_lists = output.groupby(grouping_columns)["Allele"].agg(_unique_join)
    allele_counts = output.groupby(grouping_columns)["Allele"].nunique()
    output["Allele Count"] = output.set_index(grouping_columns).index.map(allele_counts).astype("Int64")
    output["Alleles"] = output.set_index(grouping_columns).index.map(allele_lists)
    output["Classification (conservation)"] = output["Protein_ID"].map(classification_map).fillna("Unclassified")
    columns = [
        "Protein_ID",
        "Allele",
        "Epitope",
        "Allele Count",
        "Elution Score",
        "Rank",
        "Alleles",
        "Classification (conservation)",
        "start",
        "end",
    ]
    output = output[columns]
    return output.sort_values(by=["Allele Count", "Rank", "Protein_ID", "Epitope", "Allele"], ascending=[False, True, True, True, True]).reset_index(drop=True)


def build_mhcii_final_output(frame: pd.DataFrame, classification_map: dict[str, str]) -> pd.DataFrame:
    output = frame.copy()
    output["score"] = pd.to_numeric(output["score"], errors="coerce")
    output["rank"] = pd.to_numeric(output["rank"], errors="coerce")


    output["core_key"] = output["core_peptide"].fillna("").astype(str)
    output.loc[output["core_key"].str.strip() == "", "core_key"] = output["peptide"].astype(str)

    grouped_rows: list[dict[str, object]] = []
    for (protein_id, allele, peptide, core_key), group in output.groupby(
        ["protein_id", "allele", "peptide", "core_key"],
        dropna=False,
    ):
        grouped_rows.append(
            {
                "Protein_ID": protein_id,
                "Epitope": peptide,
                "Core": core_key,
                "Allele": allele,
                "Elution Score": group["score"].max(),
                "Rank": group["rank"].min(),
                "start": group["start"].dropna().iloc[0] if group["start"].notna().any() else pd.NA,
                "end": group["end"].dropna().iloc[0] if group["end"].notna().any() else pd.NA,
            }
        )

    output = pd.DataFrame(grouped_rows)
    # Counted per protein and epitope, not per core. One epitope can be reported against several cores.
    grouping_columns = ["Protein_ID", "Epitope"]
    allele_lists = output.groupby(grouping_columns)["Allele"].agg(_unique_join)
    allele_counts = output.groupby(grouping_columns)["Allele"].nunique()
    output["Allele Count"] = output.set_index(grouping_columns).index.map(allele_counts).astype("Int64")
    output["Alleles"] = output.set_index(grouping_columns).index.map(allele_lists)
    output["Classification (conservation)"] = output["Protein_ID"].map(classification_map).fillna("Unclassified")
    columns = [
        "Protein_ID",
        "Allele",
        "Epitope",
        "Core",
        "Allele Count",
        "Elution Score",
        "Rank",
        "Alleles",
        "Classification (conservation)",
        "start",
        "end",
    ]
    output = output[columns]
    if output.empty:
        return output
    return output.sort_values(by=["Allele Count", "Rank", "Protein_ID", "Core", "Allele"], ascending=[False, True, True, True, True]).reset_index(drop=True)


def build_mhci_compatibility_detail(frame: pd.DataFrame, classification_map: dict[str, str]) -> pd.DataFrame:
    detailed = frame.copy()
    detailed["Rank"] = pd.to_numeric(detailed["rank"], errors="coerce")
    detailed["Elution Score"] = pd.to_numeric(detailed["score"], errors="coerce")
    counts = detailed.groupby(["protein_id", "peptide"])["allele"].nunique()
    detailed["Allele Count"] = detailed.set_index(["protein_id", "peptide"]).index.map(counts)
    detailed["Classification (conservation)"] = detailed["protein_id"].map(classification_map).fillna("Unclassified")
    return detailed.rename(
        columns={"protein_id": "Protein_ID", "peptide": "Epitope", "allele": "Allele"}
    )[
        [
            "Protein_ID",
            "Epitope",
            "Classification (conservation)",
            "Allele",
            "Allele Count",
            "start",
            "end",
            "Elution Score",
            "Rank",
        ]
    ].sort_values(by=["Protein_ID", "Epitope", "Rank"], ascending=[True, True, True]).reset_index(drop=True)


def build_mhcii_compatibility_detail(frame: pd.DataFrame, classification_map: dict[str, str]) -> pd.DataFrame:
    detailed = frame.copy()
    detailed["Rank"] = pd.to_numeric(detailed["rank"], errors="coerce")
    detailed["Elution Score"] = pd.to_numeric(detailed["score"], errors="coerce")
    counts = detailed.groupby(["protein_id", "peptide", "core_peptide"])["allele"].nunique()
    detailed["Allele Count"] = detailed.set_index(["protein_id", "peptide", "core_peptide"]).index.map(counts)
    detailed["Classification (conservation)"] = detailed["protein_id"].map(classification_map).fillna("Unclassified")
    return detailed.rename(
        columns={"protein_id": "Protein_ID", "peptide": "Epitope", "core_peptide": "Core", "allele": "Allele"}
    )[
        [
            "Protein_ID",
            "Epitope",
            "Core",
            "Classification (conservation)",
            "Allele",
            "Allele Count",
            "start",
            "end",
            "Elution Score",
            "Rank",
        ]
    ].sort_values(by=["Protein_ID", "Core", "Rank"], ascending=[True, True, True]).reset_index(drop=True)


def build_compatibility_table(
    mhci_detail: pd.DataFrame,
    mhcii_detail: pd.DataFrame,
    protein_metadata: dict[str, dict[str, str]],
    minimum_overlap: int = 7,
) -> pd.DataFrame:

    rows: list[dict[str, object]] = []
    mhci_by_protein = {protein_id: group.copy() for protein_id, group in mhci_detail.groupby("Protein_ID")}

    for _, mhcii_row in mhcii_detail.iterrows():
        protein_id = str(mhcii_row["Protein_ID"])
        mhci_candidates = mhci_by_protein.get(protein_id)
        if mhci_candidates is None:
            continue
        core = str(mhcii_row["Core"])
        species = protein_metadata.get(protein_id, {}).get("species", "Unknown")
        for _, mhci_row in mhci_candidates.iterrows():
            epitope = str(mhci_row["Epitope"])
            overlap_size = longest_contiguous_overlap(core, epitope)
            if overlap_size < minimum_overlap:
                continue
            rows.append(
                {
                    "Species": species,
                    "Protein_ID": protein_id,
                    "MHC-II Epitope": mhcii_row["Epitope"],
                    "MHC-II Core": core,
                    "MHC-II Classification (conservation)": mhcii_row["Classification (conservation)"],
                    "MHC-II Allele": mhcii_row["Allele"],
                    "MHC-II Allele Count": mhcii_row["Allele Count"],
                    "MHC-II start": mhcii_row["start"],
                    "MHC-II end": mhcii_row["end"],
                    "MHC-II Elution Score": mhcii_row["Elution Score"],
                    "MHC-II Rank": mhcii_row["Rank"],
                    "MHC-I Epitope": epitope,
                    "MHC-I Classification (conservation)": mhci_row["Classification (conservation)"],
                    "MHC-I Allele": mhci_row["Allele"],
                    "MHC-I Allele Count": mhci_row["Allele Count"],
                    "MHC-I start": mhci_row["start"],
                    "MHC-I end": mhci_row["end"],
                    "MHC-I Elution Score": mhci_row["Elution Score"],
                    "MHC-I Rank": mhci_row["Rank"],
                    "Overlap Size": overlap_size,
                }
            )

    columns = [
        "Species",
        "Protein_ID",
        "MHC-II Epitope",
        "MHC-II Core",
        "MHC-II Classification (conservation)",
        "MHC-II Allele",
        "MHC-II Allele Count",
        "MHC-II start",
        "MHC-II end",
        "MHC-II Elution Score",
        "MHC-II Rank",
        "MHC-I Epitope",
        "MHC-I Classification (conservation)",
        "MHC-I Allele",
        "MHC-I Allele Count",
        "MHC-I start",
        "MHC-I end",
        "MHC-I Elution Score",
        "MHC-I Rank",
        "Overlap Size",
    ]
    compatibility = pd.DataFrame(rows, columns=columns)
    if compatibility.empty:
        return compatibility
    return compatibility.sort_values(
        by=["Protein_ID", "MHC-II Core", "Overlap Size", "MHC-II Rank", "MHC-I Rank"],
        ascending=[True, True, False, True, True],
    ).reset_index(drop=True)


def longest_contiguous_overlap(left: str, right: str) -> int:
    best = 0
    for left_index in range(len(left)):
        for right_index in range(len(right)):
            size = 0
            while (
                left_index + size < len(left)
                and right_index + size < len(right)
                and left[left_index + size] == right[right_index + size]
            ):
                size += 1
            if size > best:
                best = size
    return best
