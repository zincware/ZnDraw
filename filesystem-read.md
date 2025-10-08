# File Browser Feature Plan

## Overview

Add a local filesystem browser to ZnDraw that allows users to browse and load files directly from the web interface, without needing to specify files via CLI arguments.

## Feature Goals

- Browse local filesystem where the ZnDraw server is running
- Load files into ZnDraw rooms via REST API
- Security-conscious design with configurable access restrictions
- Optional feature controlled by CLI flag

## Architecture

### 1. CLI Changes

**File:** `src/zndraw/cli.py`

Add new CLI option:

```python
file_browser: bool = typer.Option(
    False,
    "--file-browser/--no-file-browser",
    help="Enable local filesystem browser endpoint",
),
file_browser_root: str = typer.Option(
    None,
    "--file-browser-root",
    help="Root directory for file browser (defaults to current working directory)",
),
```

Pass configuration to Flask app:

```python
flask_app.config["FILE_BROWSER_ENABLED"] = file_browser
flask_app.config["FILE_BROWSER_ROOT"] = file_browser_root or os.getcwd()
```

### 2. Backend REST API

**New file:** `src/zndraw/app/file_browser.py`

#### Endpoint Naming Decision

- `/api/file-browser` - matches CLI flag, clear purpose
- frontend route: `/file-browser`

#### REST Endpoints

##### 2.1 List Directory Contents

```
GET /api/file-browser/list?path=<relative-path>
```

**Query Parameters:**
- `path` (optional): Relative path from browser root. Defaults to root.

**Response:**
```json
{
  "current_path": "/path/to/directory",
  "items": [
    {
      "name": "example.xyz",
      "type": "file",
      "size": 1024,
      "modified": "2025-01-08T10:30:00Z",
      "supported": true
    },
    {
      "name": "subdirectory",
      "type": "directory",
      "size": null,
      "modified": "2025-01-08T09:15:00Z"
    }
  ],
  "parent": "/path/to"
}
```

**Security:**
- Validate path doesn't escape browser root (path traversal protection)
- Return 403 if file_browser is disabled
- Return 400 for invalid paths
- Filter hidden files (optional, configurable)

##### 2.2 Load File into ZnDraw

```
POST /api/file-browser/load
```

**Request Body:**
```json
{
  "path": "relative/path/to/file.xyz",
  "room": "optional-room-name",
  "start": null,
  "stop": null,
  "step": null,
  "make_default": false
}
```

**Response:**
```json
{
  "status": "queued",
  "room": "relative_path_to_file_xyz",
  "message": "File loading queued",
  "task_id": "celery-task-id"
}
```

**Behavior:**
- If `room` not provided, generate from filename using `path_to_room()`
- Call `read_file.delay()` with server URL from config
- Return task ID for tracking
- Return 403 if file_browser disabled
- Return 404 if file doesn't exist
- Return 400 if file type not supported

##### 2.3 Get Supported File Types

```
GET /file-browser/supported-types
```

**Response:**
```json
{
  "extensions": [".xyz", ".pdb", ".h5md", ".traj", ".h5", ".extxyz"],
  "descriptions": {
    ".xyz": "XYZ coordinate file",
    ".pdb": "Protein Data Bank file",
    ".h5md": "H5MD molecular data",
    ".traj": "Trajectory file",
    ".h5": "H5MD molecular data",
    ".extxyz": "Extended XYZ file"
  }
}
```

**Behavior:**
- Return list of file extensions supported by `ase.io` and `znh5md`
- Used by frontend to highlight loadable files

### 3. Implementation Details

#### Security Considerations

1. **Path Traversal Protection:**
   ```python
   def validate_path(requested_path: str, root: str) -> Path | None:
       """Validate path doesn't escape root directory."""
       root = Path(root).resolve()
       target = (root / requested_path).resolve()

       if not target.is_relative_to(root):
           return None
       return target
   ```

2. **File Type Validation:**
   - Only allow reading of supported molecular file formats
   - Reject symbolic links (optional)
   - Size limits for preview (optional)

3. **Access Control:**
   - Feature disabled by default
   - Requires explicit `--file-browser` flag
   - Optional authentication (future enhancement)

#### File Format Detection

```python
def is_supported_file(filepath: Path) -> bool:
    """Check if file is supported by ASE or znh5md."""
    supported_extensions = {
        '.xyz', '.pdb', '.cif', '.mol', '.sdf',
        '.h5md', '.h5', '.hdf5',
        '.traj', '.nc', '.extxyz'
    }
    return filepath.suffix.lower() in supported_extensions
```

#### Error Handling

- 403 Forbidden: Feature disabled or path outside root
- 404 Not Found: File/directory doesn't exist
- 400 Bad Request: Invalid path format or unsupported file type
- 500 Internal Server Error: File reading errors

### 4. Frontend Considerations

**New Component:** File Browser UI (implementation details TBD)

Features:
- Tree view or breadcrumb navigation
- File type icons
- File size and modification date display
- Click to load file into room
- Visual indication of supported vs unsupported files
- Loading status indicator

**Suggested Route:** `/file-browser`

### 5. Configuration Options

| Option | Default | Description |
|--------|---------|-------------|
| `--file-browser` | `False` | Enable file browser |
| `--file-browser-root` | `cwd` | Root directory for browsing |

### 6. Implementation Phases

#### Phase 1: Backend API (Minimal Viable Product)
- [ ] Add CLI flags for file browser
- [ ] Create `file_browser.py` with REST endpoints
- [ ] Implement path validation and security checks
- [ ] Implement `/file-browser/list` endpoint
- [ ] Implement `/file-browser/load` endpoint
- [ ] Add unit tests for path validation
- [ ] Add integration tests for endpoints

#### Phase 2: Enhanced Backend
- [ ] Implement `/file-browser/supported-types` endpoint
- [ ] Add file type detection and filtering
- [ ] Add pagination for large directories
- [ ] Add search/filter capability
- [ ] Add error handling and logging

#### Phase 3: Frontend Integration
- [ ] Design file browser UI component, use mui-x TreeView
- [ ] Implement directory navigation
- [ ] Add file loading functionality
- [ ] Show loading progress/status
- [ ] Add keyboard navigation


### 7. Testing Strategy

#### Unit Tests
- Path validation (traversal attacks)
- File type detection
- Room name generation
- Error handling

#### Integration Tests
- List directory endpoint
- Load file endpoint
- Disabled feature returns 403
- Invalid paths return appropriate errors

#### Security Tests
- Path traversal attempts (`../../../etc/passwd`)
- Symlink handling
- Large directory handling
- Concurrent requests


### 10. Open Questions

1. Should we cache directory listings? response: NO
2. Should we support file upload (opposite direction)? response: No
3. Should we allow creating new rooms from the file browser? response: yes, when selecting a file we must create a room to upload the file to.
4. Should we show file loading progress in the browser UI? response: yes
5. Should we support watching directories for new files? response: no
6. Should we integrate with file system change notifications? response: no

### 11. Dependencies

No new dependencies required. Uses existing:
- Flask for REST endpoints
- Celery for async file loading
- pathlib for path handling
- ase.io for file format detection
