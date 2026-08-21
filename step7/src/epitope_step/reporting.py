from __future__ import annotations

from .models import WorkflowResult


def format_summary(results: list[WorkflowResult]) -> str:
    lines = ["Execution summary:"]
    for result in results:
        artifact_list = ", ".join(path.name for path in result.artifacts)
        metric_list = ", ".join(f"{key}={value}" for key, value in result.metrics.items())
        lines.append(f"- {result.workflow}: {metric_list}")
        lines.append(f"  artifacts: {artifact_list}")
    return "\n".join(lines)
