from __future__ import annotations

from .models import ExecutionResult


def format_execution_summary(result: ExecutionResult) -> str:
    lines = [f"Workflow: {result.workflow}"]

    if result.written_files:
        lines.append("Written files:")
        for name, path in result.written_files.items():
            lines.append(f"- {name}: {path}")

    if result.row_counts:
        lines.append("Row counts:")
        for name, count in result.row_counts.items():
            lines.append(f"- {name}: {count}")

    if result.classification_counts:
        lines.append("Classification counts:")
        for name, count in result.classification_counts.items():
            lines.append(f"- {name}: {count}")

    return "\n".join(lines)
