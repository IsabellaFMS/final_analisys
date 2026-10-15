from __future__ import annotations

import pandas as pd


def annotate_dataframe(df: pd.DataFrame, cons_map: dict[str, str]) -> pd.DataFrame:
    
    if df.empty:
        return df.copy()

    out_df = df.copy()
    first_col = out_df.columns[0]
    out_df[first_col] = out_df[first_col].astype(str).str.strip()
    out_df["classification"] = out_df[first_col].map(cons_map).fillna("Unclassified")
    return out_df


def merge_analytics(signal_df: pd.DataFrame, deep_df: pd.DataFrame) -> pd.DataFrame:

    
    signal_columns = ["id", "prediction", "confidence_level"]
    deep_columns = ["id", "signal_present", "type"]

    s_safe = signal_df[signal_columns] if not signal_df.empty else pd.DataFrame(columns=signal_columns)
    d_safe = deep_df[deep_columns] if not deep_df.empty else pd.DataFrame(columns=deep_columns)
    merged = pd.merge(s_safe, d_safe, on="id", how="outer")

    def resolve_hierarchy(row: pd.Series) -> str:
        signalp_positive = str(row.get("prediction", "")) == "SP"
        deeptmhmm_positive = str(row.get("signal_present", "")) == "Yes"

        if signalp_positive and deeptmhmm_positive:
            return "Group 1: High Confidence Consensus (SP+/THMM+)"
        if signalp_positive and not deeptmhmm_positive:
            return "Group 2: SignalP Only (SP+/THMM-)"
        if not signalp_positive and deeptmhmm_positive:
            return "Group 3: DeepTMHMM Only (SP-/THMM+)"
        return "Group 4: Negative Consensus (SP-/THMM-)"

    merged["strategic_classification"] = merged.apply(resolve_hierarchy, axis=1)
    return merged
