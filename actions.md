## ZnDraw Worker and Task Management Architecture

This document outlines the architecture for a robust, scalable system for managing and executing tasks from different providers. The system is designed around a central backend that orchestrates work between default (Celery) workers and custom, remote Python (`ZnDraw`) workers.

The core principles are:

  - **Centralized State:** The backend, using Redis, is the single source of truth for all registered extensions, schemas, and worker availability.
  - **Event-Driven Logic:** The backend is entirely event-driven. It reacts to job submissions, worker connections, and task completions without any continuous polling, making it highly efficient.
  - **Decoupled Components:** The frontend, backend, and workers are decoupled. The frontend doesn't know who will run a job, and workers don't know about the backend's internal state.

### Component Overview

```mermaid
graph TD
    subgraph "User's Browser"
        Frontend[React App w/ TanStack Query]
    end

    subgraph "Server Infrastructure"
        Backend[Flask REST API & Socket.IO]
        CeleryWorker["Default Worker (Celery)"]
    end

    subgraph "Remote Machine"
        CustomWorker["Custom Worker (PythonClient)"]
    end

    subgraph "Data Store"
        Redis[(Redis)]
    end

    Frontend -- "HTTP API Calls (fetch schemas, submit actions)" --> Backend
    Backend -- "Manages State & Queues Jobs" --> Redis
    CustomWorker -- "HTTP API Calls (register)" --> Backend
    CustomWorker -- "Socket.IO Events ('worker_ready', 'task_complete')" --> Backend
    
    Backend -- "Push Notifications ('schemas_invalidated', 'queue_updated', 'run_task')" --> Frontend & CustomWorker
    Backend -- "Enqueues Jobs" --> CeleryWorker
```

**NOTE**: We must prohibit creating a room called `public`!

-----

### 1\. Worker Registration and Schema Validation

A custom worker registers its capabilities (extensions) via a REST API call. The backend is responsible for validating the extension, computing a schema hash to detect conflicts, and storing its metadata in Redis. This process supports multiple workers registering the same extension.

```mermaid
sequenceDiagram
    participant CustomWorker as "Custom Worker"
    participant Backend
    participant Redis
    participant Frontend

    Note over CustomWorker: User calls vis.register(MyExtension, public=True)

    CustomWorker->>+Backend: 1. POST /api/extensions (payload: {name, type, schema, public})
    Backend->>Backend: 2. Server computes SHA-256 hash of the received schema
    
    Note over Backend: Scope is determined by the 'public' flag (e.g., schemas:public or schemas:room123)
    Backend->>+Redis: 3. HGET schemas:{scope}:hashes, "MyExtension"
    Redis-->>Backend: returns existing hash (or null)
    
    alt Extension with same name and different schema exists
        Backend->>Backend: 4a. Compare SERVER-COMPUTED hash with stored hash
        Backend-->>CustomWorker: 409 Conflict Error: "Schema conflict for this extension"
    else New or matching schema
        Note over Backend: 4b. No conflict, proceeding.
        Backend->>+Redis: HSET schemas:{scope}, "MyExtension", '{schema...}'
        Backend->>+Redis: HSET schemas:{scope}:hashes, "MyExtension", "new_hash"
        Backend->>+Redis: SADD workers:for:{scope}:MyExtension, "worker_sid"
        Backend->>+Redis: SADD worker:tasks:{sid}, "MyExtension"
        Redis-->>Backend: OK
        Backend-->>CustomWorker: 201 Created or 200 OK
        
        Note over Backend: 5. Registration successful, now notify clients.
        Backend->>-Frontend: 6. emit('schemas_invalidated', {type: "analyses"}) via Socket.IO
    end
```

-----

### 2\. Job Lifecycle: Routing and Queuing

The job lifecycle is managed through an event-driven workflow, ensuring efficiency and providing real-time feedback to the user.

#### Job Dispatching Logic

When a job is submitted, the backend first determines the provider responsible for it (a custom worker or a default Celery worker) by checking for a registration first in the room-specific scope, then in the public scope. It then routes the task accordingly.

```mermaid
sequenceDiagram
    participant Frontend
    participant Backend
    participant Redis
    participant CustomWorker as "Custom Worker"
    participant CeleryWorker as "Default Worker"

    Note over Backend, Redis: On server startup, default jobs are pre-registered with provider="celery".

    Frontend->>+Backend: POST /api/actions (Submit 'SelectAll' job)
    Backend->>+Redis: HGET extensions:room123, "SelectAll" (checks room-specific first)
    Redis-->>Backend: returns null
    Backend->>+Redis: HGET extensions:public, "SelectAll" (falls back to public)
    Redis-->>Backend: returns '{"provider": "celery"}'
    
    alt Provider is "celery"
        Note over Backend: Routing to default worker.
        Backend->>CeleryWorker: Enqueue task in Celery's queue
    else Provider is a client SID
        Note over Backend: Routing to custom worker via queuing logic.
        Backend->>Backend: Engage Custom Worker Task Handling
    end
    
    Backend-->>Frontend: HTTP Response: {"status": "submitted"}
```

#### Custom Worker Task Handling and Queuing

If a job is routed to a custom worker, the backend uses a Redis-based queuing system to manage the workload and dispatch tasks to available workers.

```mermaid
sequenceDiagram
    participant Frontend
    participant Backend
    participant Redis
    participant CustomWorker as "Custom Worker"

    Note over CustomWorker, Backend: Worker connects and registers its extensions.
    CustomWorker->>+Backend: emit('worker_ready')
    Note over Backend: Server knows worker's capabilities from registration.
    Backend->>+Redis: SADD idle_workers:Analysis, "worker_sid"

    %% --- Job Arrives, Worker is Available ---
    Frontend->>+Backend: POST /api/actions (Submit "Job 1")
    Backend->>+Redis: SPOP idle_workers:Analysis
    Redis-->>Backend: returns "worker_sid"
    Backend->>-CustomWorker: emit('run_task', "Job 1")
    Backend-->>Frontend: HTTP Response: {"status": "submitted"}
    Backend->>-Frontend: emit('queue_updated', {task_name: 'Analysis', length: 0})

    %% --- Job Arrives, Worker is Busy ---
    Frontend->>+Backend: POST /api/actions (Submit "Job 2")
    Backend->>+Redis: SPOP idle_workers:Analysis
    Redis-->>Backend: returns null (no idle workers)
    Backend->>+Redis: LPUSH queue:task:Analysis, "Job 2"
    Redis-->>Backend: returns new queue length (e.g., 1)
    Backend-->>Frontend: HTTP Response: {"status": "queued", "position": 1}
    Backend->>-Frontend: emit('queue_updated', {task_name: 'Analysis', length: 1})
    
    %% --- Worker Finishes Job, Pulls Next from Queue ---
    CustomWorker->>+Backend: emit('task_complete', {task_name: 'Analysis'})
    Backend->>Redis: RPOP queue:task:Analysis
    Redis-->>Backend: returns "Job 2"
    Backend->>CustomWorker: emit('run_task', "Job 2")
    Backend->>Frontend: emit('queue_updated', {task_name: 'Analysis', length: 0})
```

-----

### 3\. Worker Disconnect and State Cleanup

To maintain system integrity, the backend must automatically clean up all resources associated with a worker if its connection is lost. This prevents "ghost" workers and ensures the UI accurately reflects available capabilities.

```mermaid
sequenceDiagram
    participant CustomWorker as "Custom Worker"
    participant Backend
    participant Redis
    participant Frontend

    CustomWorker--xBackend: Connection Lost
    Note over Backend: 'disconnect' event fires for "worker_sid".

    Backend->>+Redis: 1. SMEMBERS worker:tasks:worker_sid
    Redis-->>Backend: returns ["Analysis", "MyModifier"]
    
    Note over Backend: 2. For each task, remove worker from pools.
    Backend->>+Redis: SREM workers:for:Analysis, "worker_sid"
    Backend->>+Redis: SREM idle_workers:Analysis, "worker_sid"
    Backend->>+Redis: SREM workers:for:MyModifier, "worker_sid"
    Backend->>+Redis: SREM idle_workers:MyModifier, "worker_sid"
    
    Backend->>+Redis: 3. Delete the worker's own capability list.
    Redis->>Redis: DEL worker:tasks:worker_sid
    
    Note over Backend: 4. Notify all clients of the change.
    Backend->>Frontend: 5. emit('schemas_invalidated', ...)
```

-----

### 4\. Redis Key Structure Overview

The entire system is orchestrated using a structured set of keys in Redis to manage state efficiently.

```mermaid
graph TD
    subgraph "Schema & Metadata (HASH)"
        A["schemas:{scope}"] --> B("TaskName: '{schema...}'")
        C["schemas:{scope}:hashes"] --> D("TaskName: 'sha256_hash'")
    end

    subgraph "Worker & Task Pools (SET)"
        E["workers:for:{scope}:{TaskName}"] --> F1("sid_A")
        E --> F2("sid_B")
        G["idle_workers:{scope}:{TaskName}"] --> H("sid_A")
    end
    
    subgraph "Job Queues (LIST)"
        I["queue:task:{scope}:{TaskName}"] --> J1["Job 2"]
        J1 --> J2["Job 1"]
    end

    subgraph "Reverse Mapping (SET)"
        K["worker:tasks:{sid}"] --> L1("{TaskNameA}")
        K --> L2("{TaskNameB}")
    end
```

#  Diagram 1: Worker Registration & Initial State

```mermaid
sequenceDiagram
    participant Worker
    participant Server
    participant Redis
    participant UI

    Worker->>Server: Connect (Socket.IO)
    Server->>Worker: Connection Acknowledged (SID Assigned)
    Worker->>Server: register_worker(worker_id, <sio_sid>)
    Server->>Redis: HSET worker:<worker_id> status:idle, sid:<sio_sid>
    Server->>Redis: HSET sid_to_worker <sio_sid>:<worker_id>
    Server->>UI: ui_update_worker({id:worker_id, status:idle}) (Socket.IO)
    Note over Worker,UI: Worker is now idle and registered. UI shows it online.
```

# Diagram 2: Job Submission & Worker Processing (Successful)

```mermaid
sequenceDiagram
    participant Client
    participant Server
    participant Redis
    participant Worker
    participant UI

    Client->>Server: POST /api/jobs (Job Data)
    Server->>Redis: LPUSH job_queue (Serialized Job)
    Server->>Worker: new_job_available (Socket.IO - Optional Notification)
    Note over Server: Server can immediately notify workers that a new job exists.
    Worker->>Server: GET /api/jobs/next?worker_id=<worker_id> (Worker Poll)
    Server->>Redis: RPOPLPUSH job_queue processing:<worker_id>
    alt Job Available
        Redis-->>Server: Job Data
        Server-->>Worker: Job Data (HTTP 200)
        Worker->>Worker: Process Job...
        Worker->>Server: update_worker_state({worker_id, status:busy, job_id}) (Socket.IO)
        Server->>Redis: HSET worker:<worker_id> status:busy, job_id:<job_id>
        Server->>UI: ui_update_worker({id:worker_id, status:busy, job_id}) (Socket.IO)
        Worker->>Server: POST /api/jobs/result/<job_id> (Job Result)
        Server->>Redis: Remove job from processing:<worker_id> list (or mark as complete)
        Worker->>Server: update_worker_state({worker_id, status:idle}) (Socket.IO)
        Server->>Redis: HSET worker:<worker_id> status:idle, job_id:null
        Server->>UI: ui_update_worker({id:worker_id, status:idle}) (Socket.IO)
        Server-->>Client: Job Queued (HTTP 202)
    else No Job Available
        Redis-->>Server: null
        Server-->>Worker: No Content (HTTP 204)
        Worker->>Worker: Wait for next poll/notification
    end
```

# Diagram 3: Worker Disconnection & Job Re-queue (Fault Tolerance)

```mermaid
sequenceDiagram
    participant Worker
    participant Server
    participant Redis
    participant UI

    Worker--xServer: Unexpected Disconnect (e.g., Crash, Network Loss)
    Server->>Redis: HGET sid_to_worker <disconnected_sid>
    Redis-->>Server: <worker_id>
    Server->>Redis: HSET worker:<worker_id> status:offline
    Server->>Redis: HDEL sid_to_worker <disconnected_sid>
    Server->>Redis: RPOPLPUSH processing:<worker_id> job_queue
    Note over Server,Redis: Atomically moves the job back to the main queue for another worker.
    Server->>UI: ui_update_worker({id:worker_id, status:offline}) (Socket.IO)
    Note over Worker,UI: The job is now safely back in the queue, and the UI reflects the worker's status.
```

# Diagram 4: UI Updates (Simplified)

```mermaid
sequenceDiagram
    participant UI_Client
    participant Server
    participant Redis
    participant Worker

    UI_Client->>Server: Connect (Socket.IO)
    UI_Client->>Server: subscribe_ui() (Socket.IO)
    Server->>UI_Client: Connection Acknowledged (SID Assigned), joined 'ui' room
    Note over Server: (Later, when a worker state changes)
    Worker->>Server: update_worker_state(...)
    Server->>Redis: HSET worker:<worker_id> status:<new_status>
    Server->>UI_Client: ui_update_worker({id:worker_id, status:<new_status>, ...}) (Socket.IO, to 'ui' room)
    Note over UI_Client: UI updates in real-time.
```