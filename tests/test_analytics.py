"""Tests for the analytics module."""

from datetime import datetime, timedelta

import pytest

from zndraw.app.analytics import aggregate_job_stats, get_daily_stats


def test_aggregate_job_stats_empty():
    """Test aggregate_job_stats with empty job list."""
    stats = aggregate_job_stats([])

    assert stats.total_jobs == 0
    assert stats.completed == 0
    assert stats.failed == 0
    assert stats.success_rate == 0.0
    assert stats.avg_wait_time_ms == 0.0
    assert stats.avg_execution_time_ms == 0.0
    assert stats.last_used is None
    assert stats.error_breakdown == []


def test_aggregate_job_stats_single_completed_job():
    """Test aggregate_job_stats with a single completed job."""
    jobs = [
        {
            "id": "job1",
            "status": "completed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T12:00:00",
        }
    ]

    stats = aggregate_job_stats(jobs)

    assert stats.total_jobs == 1
    assert stats.completed == 1
    assert stats.failed == 0
    assert stats.success_rate == 100.0
    assert stats.avg_wait_time_ms == 1000.0
    assert stats.avg_execution_time_ms == 2000.0
    assert stats.last_used == "2024-01-01T12:00:00"
    assert stats.error_breakdown == []


def test_aggregate_job_stats_multiple_jobs():
    """Test aggregate_job_stats with multiple jobs."""
    jobs = [
        {
            "id": "job1",
            "status": "completed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T12:00:00",
        },
        {
            "id": "job2",
            "status": "completed",
            "wait_time_ms": "3000",
            "execution_time_ms": "4000",
            "completed_at": "2024-01-01T13:00:00",
        },
        {
            "id": "job3",
            "status": "failed",
            "wait_time_ms": "500",
            "execution_time_ms": "1000",
            "completed_at": "2024-01-01T14:00:00",
            "error": "Something went wrong",
        },
    ]

    stats = aggregate_job_stats(jobs)

    assert stats.total_jobs == 3
    assert stats.completed == 2
    assert stats.failed == 1
    assert stats.success_rate == pytest.approx(66.67, rel=0.01)
    assert stats.avg_wait_time_ms == pytest.approx(1500.0)
    assert stats.avg_execution_time_ms == pytest.approx(2333.33, rel=0.01)
    assert stats.last_used == "2024-01-01T14:00:00"
    assert len(stats.error_breakdown) == 1
    assert stats.error_breakdown[0]["error"] == "Something went wrong"
    assert stats.error_breakdown[0]["count"] == 1


def test_aggregate_job_stats_with_queued_jobs():
    """Test that queued jobs don't affect averages."""
    jobs = [
        {
            "id": "job1",
            "status": "completed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T12:00:00",
        },
        {
            "id": "job2",
            "status": "queued",
            "wait_time_ms": "",
            "execution_time_ms": "",
        },
    ]

    stats = aggregate_job_stats(jobs)

    assert stats.total_jobs == 2
    assert stats.completed == 1
    assert stats.failed == 0
    assert stats.success_rate == 50.0
    assert stats.avg_wait_time_ms == 1000.0
    assert stats.avg_execution_time_ms == 2000.0


def test_aggregate_job_stats_error_truncation():
    """Test that long errors are truncated."""
    long_error = "x" * 200
    jobs = [
        {
            "id": "job1",
            "status": "failed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T12:00:00",
            "error": long_error,
        }
    ]

    stats = aggregate_job_stats(jobs)

    assert len(stats.error_breakdown) == 1
    assert len(stats.error_breakdown[0]["error"]) == 100
    assert stats.error_breakdown[0]["error"] == long_error[:100]


def test_aggregate_job_stats_multiple_same_errors():
    """Test error breakdown with repeated errors."""
    jobs = [
        {
            "id": "job1",
            "status": "failed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T12:00:00",
            "error": "Error A",
        },
        {
            "id": "job2",
            "status": "failed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T13:00:00",
            "error": "Error A",
        },
        {
            "id": "job3",
            "status": "failed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "completed_at": "2024-01-01T14:00:00",
            "error": "Error B",
        },
    ]

    stats = aggregate_job_stats(jobs)

    assert len(stats.error_breakdown) == 2
    # Sorted by count descending
    assert stats.error_breakdown[0]["error"] == "Error A"
    assert stats.error_breakdown[0]["count"] == 2
    assert stats.error_breakdown[1]["error"] == "Error B"
    assert stats.error_breakdown[1]["count"] == 1


def test_get_daily_stats_empty():
    """Test get_daily_stats with empty job list."""
    daily_stats = get_daily_stats([], days=7)

    assert len(daily_stats) == 7
    for day_stats in daily_stats:
        assert "date" in day_stats
        assert day_stats["job_count"] == 0
        assert day_stats["success_rate"] == 0.0
        assert day_stats["avg_wait_ms"] == 0.0
        assert day_stats["avg_exec_ms"] == 0.0


def test_get_daily_stats_single_day():
    """Test get_daily_stats with jobs from a single day."""
    today = datetime.utcnow().strftime("%Y-%m-%d")
    jobs = [
        {
            "id": "job1",
            "status": "completed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "created_at": f"{today}T12:00:00",
            "completed_at": f"{today}T12:00:10",
        },
        {
            "id": "job2",
            "status": "completed",
            "wait_time_ms": "3000",
            "execution_time_ms": "4000",
            "created_at": f"{today}T13:00:00",
            "completed_at": f"{today}T13:00:07",
        },
    ]

    daily_stats = get_daily_stats(jobs, days=3)

    assert len(daily_stats) == 3
    # Find today's stats
    today_stats = next(s for s in daily_stats if s["date"] == today)
    assert today_stats["job_count"] == 2
    assert today_stats["success_rate"] == 100.0
    assert today_stats["avg_wait_ms"] == 2000.0
    assert today_stats["avg_exec_ms"] == 3000.0


def test_get_daily_stats_multiple_days():
    """Test get_daily_stats with jobs across multiple days."""
    today = datetime.utcnow()
    yesterday = today - timedelta(days=1)

    jobs = [
        {
            "id": "job1",
            "status": "completed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "created_at": today.strftime("%Y-%m-%dT12:00:00"),
            "completed_at": today.strftime("%Y-%m-%dT12:00:10"),
        },
        {
            "id": "job2",
            "status": "completed",
            "wait_time_ms": "3000",
            "execution_time_ms": "4000",
            "created_at": yesterday.strftime("%Y-%m-%dT13:00:00"),
            "completed_at": yesterday.strftime("%Y-%m-%dT13:00:07"),
        },
    ]

    daily_stats = get_daily_stats(jobs, days=3)

    assert len(daily_stats) == 3
    # Verify each day has the correct data
    today_str = today.strftime("%Y-%m-%d")
    yesterday_str = yesterday.strftime("%Y-%m-%d")

    today_stats = next(s for s in daily_stats if s["date"] == today_str)
    assert today_stats["job_count"] == 1

    yesterday_stats = next(s for s in daily_stats if s["date"] == yesterday_str)
    assert yesterday_stats["job_count"] == 1


def test_get_daily_stats_invalid_timestamps():
    """Test get_daily_stats handles invalid timestamps gracefully."""
    jobs = [
        {
            "id": "job1",
            "status": "completed",
            "wait_time_ms": "1000",
            "execution_time_ms": "2000",
            "created_at": "invalid-timestamp",
            "completed_at": "2024-01-01T12:00:10",
        },
    ]

    # Should not raise an exception
    daily_stats = get_daily_stats(jobs, days=7)
    assert len(daily_stats) == 7


def test_aggregate_job_stats_includes_running_jobs_in_averages():
    """Test that running jobs contribute to wait_time averages."""
    jobs = [
        {
            "id": "job1",
            "status": "running",
            "wait_time_ms": "1000",
            "execution_time_ms": "",
            "created_at": "2024-01-01T12:00:00",
        },
        {
            "id": "job2",
            "status": "completed",
            "wait_time_ms": "2000",
            "execution_time_ms": "3000",
            "completed_at": "2024-01-01T13:00:00",
        },
    ]

    stats = aggregate_job_stats(jobs)

    # Both jobs should contribute to wait time average
    assert stats.total_jobs == 2
    assert stats.completed == 1
    assert stats.avg_wait_time_ms == pytest.approx(1500.0)  # (1000 + 2000) / 2
    # Only completed job contributes to execution time
    assert stats.avg_execution_time_ms == pytest.approx(3000.0)
    assert stats.success_rate == 50.0


def test_aggregate_job_stats_running_jobs_no_execution_time():
    """Test that running jobs don't affect execution time averages."""
    jobs = [
        {
            "id": "job1",
            "status": "running",
            "wait_time_ms": "500",
            "execution_time_ms": "",
            "created_at": "2024-01-01T12:00:00",
        },
        {
            "id": "job2",
            "status": "running",
            "wait_time_ms": "600",
            "execution_time_ms": "",
            "created_at": "2024-01-01T12:01:00",
        },
    ]

    stats = aggregate_job_stats(jobs)

    assert stats.total_jobs == 2
    assert stats.completed == 0
    assert stats.failed == 0
    assert stats.success_rate == 0.0
    assert stats.avg_wait_time_ms == pytest.approx(550.0)
    assert stats.avg_execution_time_ms == 0.0  # No completed jobs
