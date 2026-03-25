# tests/test_stress.py
"""Async stress tests for concurrent task operations."""

import asyncio

import pytest

from zndraw_joblib.schemas import PaginatedResponse, TaskResponse


@pytest.mark.asyncio
async def test_concurrent_task_submissions(async_client):
    """Submit 100 tasks concurrently."""
    # First create a worker and register a job
    worker_resp = await async_client.post("/v1/joblib/workers")
    assert worker_resp.status_code == 201
    worker_id = worker_resp.json()["id"]

    job_resp = await async_client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "StressJob", "worker_id": worker_id},
    )
    assert job_resp.status_code == 201

    async def submit_task(i: int):
        return await async_client.post(
            "/v1/joblib/rooms/stress_room/tasks/@global:modifiers:StressJob",
            json={"payload": {"index": i}},
        )

    # Submit 100 tasks concurrently
    responses = await asyncio.gather(*[submit_task(i) for i in range(100)])
    assert all(r.status_code == 202 for r in responses)

    # Verify all tasks exist
    list_resp = await async_client.get("/v1/joblib/rooms/stress_room/tasks?limit=200")
    assert list_resp.status_code == 200
    page = PaginatedResponse[TaskResponse].model_validate(list_resp.json())
    assert page.total == 100
    assert len(page.items) == 100


@pytest.mark.asyncio
async def test_concurrent_claims_multiple_workers(async_client):
    """10 workers claim concurrently from 50 tasks."""
    # Create 10 workers
    worker_ids = []
    for i in range(10):
        worker_resp = await async_client.post("/v1/joblib/workers")
        assert worker_resp.status_code == 201
        worker_id = worker_resp.json()["id"]
        worker_ids.append(worker_id)

        # Register same job for each worker
        await async_client.put(
            "/v1/joblib/rooms/@global/jobs",
            json={
                "category": "modifiers",
                "name": "ConcurrentJob",
                "worker_id": worker_id,
            },
        )

    # Submit 50 tasks
    for i in range(50):
        resp = await async_client.post(
            "/v1/joblib/rooms/concurrent_room/tasks/@global:modifiers:ConcurrentJob",
            json={"payload": {"index": i}},
        )
        assert resp.status_code == 202

    # All workers claim concurrently (multiple rounds)
    claimed_task_ids = set()

    async def claim_tasks(worker_id: str):
        """Claim until no more tasks."""
        local_claimed = []
        for _ in range(10):  # Each worker tries to claim up to 10 tasks
            resp = await async_client.post(
                "/v1/joblib/tasks/claim",
                json={"worker_id": worker_id},
            )
            assert resp.status_code == 200
            data = resp.json()
            if data["task"] is None:
                break
            local_claimed.append(data["task"]["id"])
        return local_claimed

    results = await asyncio.gather(*[claim_tasks(w) for w in worker_ids])

    # Collect all claimed task IDs
    for claimed_list in results:
        for task_id in claimed_list:
            claimed_task_ids.add(task_id)

    # All 50 tasks should have been claimed, no duplicates
    assert len(claimed_task_ids) == 50


@pytest.mark.asyncio
async def test_no_double_claim_under_contention(async_client):
    """Multiple workers trying to claim same task simultaneously - only 1 should succeed."""
    # Create 10 workers
    worker_ids = []
    for i in range(10):
        worker_resp = await async_client.post("/v1/joblib/workers")
        worker_id = worker_resp.json()["id"]
        worker_ids.append(worker_id)

        await async_client.put(
            "/v1/joblib/rooms/@global/jobs",
            json={
                "category": "modifiers",
                "name": "SingleJob",
                "worker_id": worker_id,
            },
        )

    # Submit only 1 task
    await async_client.post(
        "/v1/joblib/rooms/single_room/tasks/@global:modifiers:SingleJob",
        json={"payload": {"data": "unique"}},
    )

    # All 10 workers try to claim at once
    async def try_claim(worker_id: str):
        resp = await async_client.post(
            "/v1/joblib/tasks/claim",
            json={"worker_id": worker_id},
        )
        return resp.json()

    results = await asyncio.gather(*[try_claim(w) for w in worker_ids])

    # Only 1 should have claimed the task
    claimed_count = sum(1 for r in results if r["task"] is not None)
    assert claimed_count == 1


@pytest.mark.asyncio
async def test_large_queue_concurrent_submit_and_claim(async_client):
    """Submit and claim happening concurrently."""
    # Create workers first
    worker_ids = []
    for i in range(5):
        worker_resp = await async_client.post("/v1/joblib/workers")
        worker_id = worker_resp.json()["id"]
        worker_ids.append(worker_id)

        await async_client.put(
            "/v1/joblib/rooms/@global/jobs",
            json={
                "category": "modifiers",
                "name": "MixedJob",
                "worker_id": worker_id,
            },
        )

    # Track submissions and claims
    submitted = []
    claimed = []

    async def submit_batch(start: int, count: int):
        """Submit a batch of tasks."""
        for i in range(count):
            resp = await async_client.post(
                "/v1/joblib/rooms/mixed_room/tasks/@global:modifiers:MixedJob",
                json={"payload": {"batch": start, "index": i}},
            )
            if resp.status_code == 202:
                submitted.append(resp.json()["id"])

    async def claim_batch(worker_id: str, attempts: int):
        """Try to claim multiple tasks."""
        for _ in range(attempts):
            resp = await async_client.post(
                "/v1/joblib/tasks/claim",
                json={"worker_id": worker_id},
            )
            if resp.status_code == 200 and resp.json()["task"] is not None:
                claimed.append(resp.json()["task"]["id"])
            await asyncio.sleep(0.01)  # Small delay to allow submitters to catch up

    # Run submitters and claimers concurrently
    tasks = [
        submit_batch(0, 30),
        submit_batch(100, 30),
        claim_batch(worker_ids[0], 20),
        claim_batch(worker_ids[1], 20),
        claim_batch(worker_ids[2], 20),
        claim_batch(worker_ids[3], 20),
        claim_batch(worker_ids[4], 20),
    ]

    await asyncio.gather(*tasks)

    # Verify no duplicate claims
    unique_claimed = set(claimed)
    assert len(unique_claimed) == len(claimed), "Some tasks were double-claimed!"

    # All claimed tasks should have been from submitted tasks
    for task_id in unique_claimed:
        assert task_id in submitted


@pytest.mark.asyncio
async def test_rapid_fire_claims_same_worker(async_client):
    """Single worker makes rapid claim requests."""
    # Create worker and register job
    worker_resp = await async_client.post("/v1/joblib/workers")
    worker_id = worker_resp.json()["id"]

    await async_client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "RapidJob", "worker_id": worker_id},
    )

    # Submit 20 tasks
    for i in range(20):
        await async_client.post(
            "/v1/joblib/rooms/rapid_room/tasks/@global:modifiers:RapidJob",
            json={"payload": {"index": i}},
        )

    # Single worker makes 30 concurrent claim requests (more than available tasks)
    async def rapid_claim():
        resp = await async_client.post(
            "/v1/joblib/tasks/claim",
            json={"worker_id": worker_id},
        )
        return resp.json()

    results = await asyncio.gather(*[rapid_claim() for _ in range(30)])

    # Count successful claims
    successful = [r for r in results if r["task"] is not None]
    failed = [r for r in results if r["task"] is None]

    # Should have claimed exactly 20 tasks (all available)
    assert len(successful) == 20
    assert len(failed) == 10  # 30 attempts - 20 tasks = 10 nulls

    # No duplicates
    claimed_ids = [r["task"]["id"] for r in successful]
    assert len(set(claimed_ids)) == 20
