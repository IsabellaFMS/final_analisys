from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def export_csv(df: pd.DataFrame, out_path: Path) -> None:
    if df.empty:
        logger.warning(f"Writing an empty dataframe to {out_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        df.to_csv(out_path, index=False, sep=";", encoding="utf-8")
        logger.info(f"Wrote output table to {out_path}")
    except Exception as exc:
        logger.error(f"Failed to write {out_path}: {exc}")
