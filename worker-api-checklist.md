Worker API Refactoring Plan

 Overview

 Refactor from REST polling to socket-based push system for job assignment. Workers will receive jobs via socket events instead of polling
 /api/jobs/next.

 State Changes

 Job States (NEW)

 - pending - Waiting for idle worker
 - assigned - Emitted to worker, awaiting confirmation
 - running - Worker actively processing
 - completed - Successfully finished
 - failed - Error or worker disconnect

 Worker States (per extension)

 - idle - Available for work
 - assigned - Job emitted, awaiting worker response (with timestamp)
 - running - Worker processing job

 Backend Changes

 1. Remove Old Polling System

 File: src/zndraw/app/job_routes.py
 - ❌ Delete get_next_job_agnostic() endpoint (lines 228-473)
 - ❌ Delete helper function _get_next_job_from_queue() if exists
 - ✅ Keep update_job_status() but modify for new flow

 2. Add Pending Job Management

 File: src/zndraw/app/redis_keys.py
 - ✅ Add pending_jobs sorted set key to ExtensionKeys (score = timestamp)
 - ✅ Add worker state hash: {base}:worker_states (worker_id → {state, assigned_at, job_id})

 3. Refactor Job Submission

 File: src/zndraw/app/extension_routes.py
 - ✅ Modify _submit_extension_impl() (lines 24-231):
   - Create job with status pending
   - Check for idle workers using SRANDMEMBER {idle_workers} 1
   - If idle worker exists:
       - Move worker: idle → assigned with HSET worker_states {worker_id} {state:assigned, assigned_at:timestamp, job_id:...}
     - Update job status to assigned
     - Emit job:assigned event to specific worker with job data
     - Remove from pending set if it was there
   - If no idle worker:
       - Add to ZADD pending_jobs {timestamp} {job_id}
   - For Celery extensions:
       - Trigger celery_job_worker.delay() immediately
     - Job starts as assigned to celery worker
   - Emit job state update to room clients
   - Don't call emit_queue_update() (removing this)

 4. Add Job Assignment Socket Event

 File: src/zndraw/app/constants.py
 - ✅ Add new event constant: JOB_ASSIGNED = "job:assigned"
 - ✅ Add JOB_STATE_CHANGED = "job:state_changed"

 File: src/zndraw/app/events.py (or new file job_dispatcher.py)
 - ✅ Add function emit_job_to_worker(worker_id, job_data)
 - ✅ Add function emit_job_state_change(room_id, job_id, status, metadata)

 5. Update Job Status Endpoint

 File: src/zndraw/app/job_routes.py
 - ✅ Modify update_job_status() (lines 136-225):
   - Accept new status: processing (in addition to completed/failed)
   - Status: processing
       - Validate job is in assigned state
     - Update job: assigned → running
     - Update worker state: assigned → running
     - Emit job state change to room
   - Status: completed/failed
       - Update job state
     - Update worker: assigned/running → idle
     - Remove from worker_states hash
     - Trigger pending job assignment:
           - Call assign_pending_jobs(room_id, category, extension, worker_id)
     - Emit job state change to room
     - Don't call emit_queue_update() (removing)

 6. Assign Pending Jobs

 File: src/zndraw/app/job_dispatcher.py (NEW)
 - ✅ Create assign_pending_jobs(room_id, category, extension, worker_id=None):
   - Get all idle workers for extension (or specific worker_id)
   - Get pending jobs from sorted set (oldest first)
   - For each idle worker:
       - Pop oldest pending job
     - Move worker to assigned state
     - Update job to assigned state
     - Emit job:assigned to worker
     - Emit job state change to room
   - Return number of assignments made

 7. Worker Registration Changes

 File: src/zndraw/app/worker_routes.py
 - ✅ Modify register_worker() (lines 21-254):
   - After adding to idle_workers
   - Call assign_pending_jobs() for each registered extension
   - Emit schema invalidation (keep existing)

 8. Worker Disconnect Handling

 File: src/zndraw/app/events.py
 - ✅ Modify handle_disconnect() (lines 251-491):
   - Check if worker has assigned/running jobs in worker_states
   - For each active job:
       - Mark job as failed with error "Worker disconnected"
     - Emit job state change to room
   - Clean up worker from all sets and states
   - Call assign_pending_jobs() for affected extensions (if other workers idle)

 9. Remove Queue Update System

 File: src/zndraw/app/queue_manager.py
 - ❌ Delete emit_queue_update() function (lines 13-77)
 - ⚠️ Remove all calls to emit_queue_update() throughout codebase

 10. Update Job Manager

 File: src/zndraw/app/job_manager.py
 - ✅ Update JobStatus enum to include PENDING, ASSIGNED
 - ✅ Modify create_job() to accept initial status (default: pending)
 - ✅ Add assign_job(job_id, worker_id) method
 - ✅ Add start_job_from_assigned() method (assigned → running transition)
 - ✅ Update all state transition methods to emit job state changes

 11. Update Celery Worker

 File: src/zndraw/app/tasks.py
 - ✅ Modify celery_job_worker() (lines 367-466):
   - Remove call to /api/jobs/next
   - Job data passed as task parameter (from .delay() call)
   - After getting job data, immediately call PUT status endpoint with processing
   - Execute extension
   - Report completion/failure (keep existing logic)

 12. Add Monitoring Endpoints (for frontend)

 File: src/zndraw/app/extension_routes.py or job_routes.py
 - ✅ Add GET /api/rooms/{room_id}/extensions/{scope}/{category}/{extension}/stats:
   - Return: {idle_workers: N, busy_workers: N, pending_jobs: N}
 - ✅ Modify job list endpoint to include queue position for pending jobs

 Python Client Changes

 13. Remove Polling, Add Socket Listener

 File: src/zndraw/api_manager.py
 - ❌ Delete get_next_job() method (lines 158-189)
 - ✅ Add socket event handler registration in register_extension():
   - Listen for job:assigned event
   - Handler validates job data and calls extension's run() method
   - Send PUT /status with processing before starting
   - Send PUT /status with completed/failed after finishing

 File: src/zndraw/extensions/abc.py
 - ✅ Modify Extension class to handle socket-based job assignment
 - ✅ Remove polling loop from extension execution

 Frontend Changes

 14. Update Socket Event Handlers

 File: app/src/hooks/useSocketManager.ts
 - ❌ Remove onQueueUpdate handler (lines 163-196)
 - ✅ Add onJobStateChanged handler:
   - Update job state in cache
   - Invalidate jobs query
 - ✅ Keep onSchemaInvalidate (worker count changes)

 15. Add Queue Position Display

 File: app/src/components/* (sidebar/job components)
 - ✅ Fetch queue position for pending jobs
 - ✅ Display: "3 jobs ahead in queue" or "Assigned to worker" or "Processing..."
 - ✅ Add progress bar with segments: pending → assigned → running → finished/failed
   - Gray for pending, yellow for assigned, blue for running, green/red for finished

 16. Add Job History

 File: app/src/components/*
 - ✅ Display recent job history per extension (clears on new submission)
 - ✅ Show states: pending, assigned (with duration), running (with duration), completed/failed

 Migration Strategy

 1. Phase 1: Backend infrastructure
   - Add new Redis keys and state tracking
   - Add socket events and job dispatcher
   - Keep old endpoints working
 2. Phase 2: Backend job flow
   - Update submit logic to emit to workers
   - Update status endpoint for new states
   - Update disconnect handling
 3. Phase 3: Remove old system
   - Remove /api/jobs/next endpoint
   - Remove queue_update emits
   - Update Celery workers
 4. Phase 4: Client updates
   - Update Python client to socket-based
   - Update frontend components
 5. Phase 5: Testing
   - Test worker registration → job assignment
   - Test pending jobs → worker becomes available
   - Test worker disconnect during assigned/running
   - Test Celery vs remote worker parity
   - Test multiple workers, multiple jobs

 Key Files to Modify

 Backend (10 files):
 1. src/zndraw/app/job_routes.py - Remove polling, update status endpoint
 2. src/zndraw/app/extension_routes.py - New submit logic with worker check
 3. src/zndraw/app/job_dispatcher.py - NEW - Pending job assignment
 4. src/zndraw/app/job_manager.py - New states and transitions
 5. src/zndraw/app/worker_routes.py - Trigger pending assignment on register
 6. src/zndraw/app/events.py - Disconnect handling, fail jobs
 7. src/zndraw/app/redis_keys.py - New keys for pending jobs, worker states
 8. src/zndraw/app/constants.py - New socket event names
 9. src/zndraw/app/queue_manager.py - DELETE queue_update
 10. src/zndraw/app/tasks.py - Update Celery worker flow

 Python Client (2 files):
 11. src/zndraw/api_manager.py - Remove polling, add socket handler
 12. src/zndraw/extensions/abc.py - Update execution model

 Frontend (2+ files):
 13. app/src/hooks/useSocketManager.ts - New event handlers
 14. app/src/components/* - Queue position, job history, progress bar

 Testing Checklist

 - Worker registers → pending job assigned immediately
 - Job submitted, no worker → job stays pending
 - Worker finishes job → next pending job assigned
 - Worker disconnects while assigned → job marked failed
 - Worker disconnects while running → job marked failed
 - Multiple workers, multiple pending jobs → correct distribution
 - Celery job flow matches remote worker flow
 - Frontend shows correct queue position
 - Frontend progress bar updates through states
 - Job history displays correctly