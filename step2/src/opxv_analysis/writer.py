from pathlib import Path

import pandas as pd


def ensure_dirs(*paths: Path) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def write_table(df: pd.DataFrame, path: Path, sep: str = "\t") -> Path:
    ensure_dirs(path.parent)
    df.to_csv(path, index=False, sep=sep)
    return path


def write_text(text: str, path: Path) -> Path:
    ensure_dirs(path.parent)
    path.write_text(text, encoding="utf-8")
    return path
