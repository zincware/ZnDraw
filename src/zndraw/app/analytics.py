"""Analytics aggregation from job data."""

import typing as t
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta

from dateutil.parser import isoparse


@dataclass
class ExtensionStats:
    """Aggregated statistics for an extension."""

    total_jobs: int
    completed: int
    failed: int
    success_rate: float
    avg_wait_time_ms: float
    avg_execution_time_ms: float
    last_used: str | None
    error_breakdown: list[dict[str, t.Any]]


def aggregate_job_stats(jobs: list[dict]) -> ExtensionStats:
    """Aggregate statistics from a list of job dictionaries.

    Args:
        jobs: List of job data dicts from Redis

    Returns:
        ExtensionStats with aggregated metrics
    """
    if not jobs:
        return ExtensionStats(
            total_jobs=0,
            completed=0,
            failed=0,
            success_rate=0.0,
            avg_wait_time_ms=0.0,
            avg_execution_time_ms=0.0,
            last_used=None,
            error_breakdown=[],
        )

    completed = sum(1 for j in jobs if j.get("status") == "completed")
    failed = sum(1 for j in jobs if j.get("status") == "failed")

    # Calculate averages (for any job with valid times, including running jobs)
    wait_times = []
    for j in jobs:
        if j.get("wait_time_ms"):
            try:
                wait_times.append(int(j["wait_time_ms"]))
            except (ValueError, TypeError):
                pass

    exec_times = []
    for j in jobs:
        if j.get("execution_time_ms"):
            try:
                exec_times.append(int(j["execution_time_ms"]))
            except (ValueError, TypeError):
                pass

    avg_wait = sum(wait_times) / len(wait_times) if wait_times else 0.0
    avg_exec = sum(exec_times) / len(exec_times) if exec_times else 0.0

    # Find most recent job
    completed_jobs = [j for j in jobs if j.get("completed_at")]
    last_used = None
    if completed_jobs:
        last_used = max(completed_jobs, key=lambda j: j.get("completed_at", ""))[
            "completed_at"
        ]

    # Error breakdown
    error_counts: dict[str, int] = {}
    for job in jobs:
        if job.get("status") == "failed" and job.get("error"):
            error = job["error"][:100]  # Truncate long errors
            error_counts[error] = error_counts.get(error, 0) + 1

    error_breakdown = [
        {"error": error, "count": count}
        for error, count in sorted(error_counts.items(), key=lambda x: -x[1])[:10]
    ]

    return ExtensionStats(
        total_jobs=len(jobs),
        completed=completed,
        failed=failed,
        success_rate=(completed / len(jobs) * 100) if jobs else 0.0,
        avg_wait_time_ms=avg_wait,
        avg_execution_time_ms=avg_exec,
        last_used=last_used,
        error_breakdown=error_breakdown,
    )


def get_daily_stats(jobs: list[dict], days: int) -> list[dict]:
    """Group jobs by day and compute daily statistics.

    Args:
        jobs: List of job data dicts from Redis
        days: Number of days to look back

    Returns:
        List of daily statistics, oldest to newest
    """
    daily: defaultdict[str, list] = defaultdict(list)

    for job in jobs:
        if job.get("created_at"):
            try:
                date_str = isoparse(job["created_at"]).strftime("%Y-%m-%d")
                daily[date_str].append(job)
            except Exception:
                pass

    # Generate last N days
    result = []
    for i in range(days):
        date = (datetime.utcnow() - timedelta(days=i)).strftime("%Y-%m-%d")
        day_jobs = daily.get(date, [])
        stats = aggregate_job_stats(day_jobs)

        result.append(
            {
                "date": date,
                "job_count": stats.total_jobs,
                "success_rate": stats.success_rate,
                "avg_wait_ms": stats.avg_wait_time_ms,
                "avg_exec_ms": stats.avg_execution_time_ms,
            }
        )

    return list(reversed(result))  # Oldest to newest
