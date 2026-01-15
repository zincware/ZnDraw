Developer Guide
===============

This guide provides an in-depth overview of ZnDraw's architecture, focusing on how communication works between the server, Python client, and web client. Understanding these patterns is essential for contributing to ZnDraw or building extensions.

Architecture Overview
---------------------

ZnDraw implements a sophisticated real-time communication system that enables collaborative visualization and manipulation of molecular structures. The system consists of three main components that communicate through WebSocket connections and a shared data layer.

.. mermaid::

    flowchart TB
        subgraph "Python Environment"
            PC[Python Client<br/>zndraw.ZnDraw]
            AS[ASE Structures<br/>Molecular Data]
            MOD[Modifiers<br/>Analysis Tasks]
        end

        subgraph "Server Layer"
            FS[Flask-SocketIO Server<br/>zndraw.app]
            CEL[Celery Workers<br/>Background Tasks]
            ZS[ZnSocket Storage<br/>Data Synchronization]
        end

        subgraph "Web Environment"
            WC[Web Client<br/>React + Three.js]
            UI[User Interface<br/>3D Visualization]
        end

        subgraph "Storage Backend"
            RED[Redis/Storage<br/>Persistent Data]
        end

        PC <--> |Socket.IO<br/>WebSocket| FS
        WC <--> |Socket.IO<br/>WebSocket| FS
        PC <--> |ZnSocket Client<br/>Data Sync| ZS
        WC <--> |ZnSocket Client<br/>Data Sync| ZS
        FS <--> |Storage Layer| ZS
        ZS <--> |Persistence| RED
        FS <--> |Task Queue| CEL

        AS --> PC
        PC --> MOD
        MOD --> PC

Core Communication Patterns
---------------------------

ZnSocket Communication Layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ZnSocket provides the foundation for real-time data synchronization. It implements a Redis-compatible API that automatically synchronizes data structures between Python and JavaScript clients.

.. mermaid::

    sequenceDiagram
        participant PC as Python Client
        participant ZS as ZnSocket Server
        participant WC as Web Client
        participant ST as Storage

        Note over PC, WC: Initial Connection Setup
        PC->>ZS: Connect znsocket.Client
        WC->>ZS: Connect JavaScript client
        ZS->>ST: Initialize storage backend

        Note over PC, WC: Data Structure Creation
        PC->>ZS: Create znsocket.List("room:token:frames")
        ZS->>ST: Store list metadata
        ZS-->>WC: Notify structure available

        Note over PC, WC: Real-time Synchronization
        PC->>ZS: Append frame data
        ZS->>ST: Update persistent storage
        ZS-->>WC: Trigger refresh callback
        WC->>WC: Update 3D visualization

        Note over PC, WC: Bidirectional Updates
        WC->>ZS: Update selection data
        ZS->>ST: Store selection state
        ZS-->>PC: Trigger refresh callback
        PC->>PC: Process modifier queue

Room-Based Data Isolation
~~~~~~~~~~~~~~~~~~~~~~~~~

ZnDraw uses a token-based room system to isolate data between different sessions. Each room maintains its own set of data structures.

.. mermaid::

    flowchart LR
        subgraph "Default Room"
            DF[room:default:frames]
            DC[room:default:config]
            DS[room:default:selection]
        end

        subgraph "Private Room A"
            AF[room:token-a:frames]
            AC[room:token-a:config]
            AS[room:token-a:selection]
        end

        subgraph "Private Room B"
            BF[room:token-b:frames]
            BC[room:token-b:config]
            BS[room:token-b:selection]
        end

        DF -.->|copy on first access| AF
        DC -.->|copy on first access| AC
        DS -.->|copy on first access| AS

        DF -.->|copy on first access| BF
        DC -.->|copy on first access| BC
        DS -.->|copy on first access| BS

Data Structure Organization
---------------------------

Frame Storage with znsocket.List
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Molecular structures are stored as frames in a ``znsocket.List``, where each frame is a ``znsocket.Dict`` containing all atomic information.

.. mermaid::

    flowchart TD
        %% Frame List
        subgraph RoomFrames["znsocket.List: room:token:frames"]
            F0["Frame 0<br>znsocket.Dict"]
            F1["Frame 1<br>znsocket.Dict"]
            F2["Frame 2<br>znsocket.Dict"]
            FN["Frame N<br>znsocket.Dict"]
        end

        %% Frame Structure
        subgraph FrameStruct["Frame Structure"]
            direction TB
            POS["positions: Array of [x,y,z]"]
            NUM["numbers: Array of ints"]
            CEL["cell: Array of [a,b,c]"]
            CON["connectivity: Array of [i,j]"]
            ARR["arrays: Dict of custom props"]
        end

        %% Connections
        F0 --> POS
        F0 --> NUM
        F0 --> CEL
        F0 --> CON
        F0 --> ARR


Configuration Management
~~~~~~~~~~~~~~~~~~~~~~~~

The configuration system uses nested ``znsocket.Dict`` structures to organize different categories of settings.

.. mermaid::

    flowchart LR
        subgraph "znsocket.Dict: room:token:config"
            direction TB
            ROOT[Config Root]

            subgraph "Particle Settings"
                PP[particle_size: Float]
                PC[particle_color: String]
                PR[particle_radius: Float]
            end

            subgraph "Bond Settings"
                BP[bond_style: String]
                BC[bond_color: String]
                BR[bond_radius: Float]
            end

            subgraph "Camera Settings"
                CP[camera_position: Array]
                CT[camera_target: Array]
                CU[camera_up: Array]
            end

            ROOT --> PP
            ROOT --> PC
            ROOT --> PR
            ROOT --> BP
            ROOT --> BC
            ROOT --> BR
            ROOT --> CP
            ROOT --> CT
            ROOT --> CU
        end

Real-time Synchronization Workflow
----------------------------------

Python to Web Client Flow
~~~~~~~~~~~~~~~~~~~~~~~~~

When the Python client updates molecular data, the changes propagate automatically to all connected web clients.

.. mermaid::

    sequenceDiagram
        participant PY as Python Client
        participant ZL as znsocket.List
        participant WC as Web Client
        participant UI as 3D Visualization

        Note over PY, UI: Adding New Molecular Structures
        PY->>PY: Load ASE structures
        PY->>ZL: vis.extend(structures)
        ZL->>ZL: Convert to znsocket.Dict frames
        ZL-->>WC: Trigger onRefresh callback
        WC->>WC: Fetch updated frame data
        WC->>UI: Update Three.js scene
        UI->>UI: Render new molecules

        Note over PY, UI: Modifier Execution
        PY->>ZL: Register modifier function
        WC->>ZL: Trigger modifier via UI
        ZL-->>PY: Modifier queue update
        PY->>PY: Execute modifier logic
        PY->>ZL: Update frames with results
        ZL-->>WC: Propagate changes
        WC->>UI: Re-render visualization

Web to Python Client Flow
~~~~~~~~~~~~~~~~~~~~~~~~~

User interactions in the web interface can trigger actions in the Python client, such as running analysis tasks.

.. mermaid::

    sequenceDiagram
        participant UI as Web UI
        participant WC as Web Client
        participant ZS as ZnSocket
        participant PY as Python Client
        participant CEL as Celery Worker

        Note over UI, CEL: User Selection and Analysis
        UI->>WC: User selects atoms
        WC->>ZS: Update selection dict
        ZS-->>PY: Selection change callback
        PY->>PY: Check modifier queue

        Note over UI, CEL: Background Task Execution
        UI->>WC: Trigger analysis task
        WC->>ZS: Add task to queue
        ZS-->>PY: Queue update callback
        PY->>CEL: Submit Celery task
        CEL->>CEL: Execute computation
        CEL->>ZS: Store results
        ZS-->>WC: Results available
        WC->>UI: Display analysis results

Callback System Implementation
------------------------------

Web Client Refresh Handlers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The web client sets up refresh callbacks for each data type to maintain real-time synchronization.

.. code-block:: typescript

    // Example: Frame synchronization
    const frames = new znsocket.List({
        client: client,
        key: `room:${token}:frames`,
    });

    frames.onRefresh(async () => {
        console.log("Frames updated externally");
        const frameCount = await frames.len();
        setTotalFrames(frameCount);
        // Trigger 3D scene update
        updateVisualization();
    });

Python Client Queue Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python client runs a background thread that continuously monitors for changes and processes modifier queues.

.. code-block:: python

    def check_queue(vis: "ZnDraw") -> None:
        """Background thread for processing modifier and public queues."""
        while True:
            process_modifier_queue(vis)
            process_public_queue(vis)
            vis.socket.sleep(1)  # Rate limiting

    # Start background processing
    vis.socket.start_background_task(check_queue, vis)
