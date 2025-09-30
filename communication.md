# Styles of communication

1. REST + SOCKET.IO
for data heavy operations
```mermaid
sequenceDiagram
    participant ClientA
    participant Server
    participant ClientB

    ClientA->>Server: REST | modify positions
    Server-->>ClientA: REST | success

    Server-->>ClientB: Socket.io | invalidate cache

    ClientB->>Server: REST | get updated positions
    Server-->>ClientB: REST | return new positions
```
2. SOCKET.IO only
for light operations
```mermaid
sequenceDiagram
    participant ClientA
    participant Server
    participant ClientB

    ClientA->>Server: Socket.io | update positions
    Server-->>ClientA: Socket.io | ack update
    Server-->>ClientB: Socket.io | broadcast new positions
    ClientB->>ClientB: Apply update directly (no REST refetch)
```

# Usage
Everything that comes from the `ase.Atoms` uses REST + SOCKET.IO.
Settings and Extensions use REST + SOCKET.IO.

currentStep uses SOCKET.IO only.
Drawing use SOCKET.IO only.
Camera use SOCKET.IO only.
Chat use SOCKET.IO + for fetching history it uses REST + SOCKET.IO.

Extension and schema management use REST + SOCKET.IO mixture (see `actions.md`).

# Nomenclature
`zndraw.extensions.analysis` or `zndraw.extensions.modifiers`
with `ExtensionType.ANALYSIS` and `ExtensionType.MODIFIER`

- Use plural for module/package names (.modifiers, .analysis).
- Use singular for enum members (ExtensionType.MODIFIER, ExtensionType.ANALYSIS).
- Use singular for extension classes (Duplicate, SelectAll).
- `category` -> Selection, Modify, Analysis, ...
- `extension` -> specific extension, e.g. "SelectAll", "Move", "Distance", ...
