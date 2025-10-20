# ZnDraw Drag & Drop / Copy-Paste / Download Feature Plan

## Overview

This document outlines a comprehensive plan for adding drag & drop file upload, copy/paste support, and file download capabilities to ZnDraw. The implementation will support both file content upload and file path detection for local file systems.

## Current State Analysis

### Existing Infrastructure
- ✅ **Celery Task System**: `read_file` task in `tasks.py` handles file loading
- ✅ **File Browser API**: `/api/file-browser/load` endpoint for local file loading
- ✅ **Format Detection**: `FORMAT_BACKENDS` maps 70+ file formats to handlers (ASE, ZnH5MD, ASE-DB)
- ✅ **Duplicate Detection**: Checks file path, size, and mtime to avoid re-uploading
- ✅ **File Upload Pattern**: Screenshot upload in `routes.py` shows multipart/form-data handling
- ✅ **Progress Tracking**: Celery task IDs can be used for progress monitoring

### Missing Features
- ❌ Drag & drop file upload
- ❌ Copy/paste file support
- ❌ File download capability
- ❌ Browser-based file upload (non-local filesystem)

---

## Feature 1: Drag & Drop Upload

### Use Cases

1. **Local Files**: User drags a file from desktop/file manager
   - Browser provides: File object with content + name
   - File content is uploaded to server
   - Server processes via Celery task

### Frontend Implementation

#### A. Canvas Drop Zone (`Canvas.tsx`)

```typescript
// Add to Canvas component wrapper
<div
  onDragOver={handleDragOver}
  onDragLeave={handleDragLeave}
  onDrop={handleDrop}
  style={{ position: 'relative', width: '100%', height: '100%' }}
>
  {isDragging && <DropOverlay />}
  <Canvas>
    {/* existing canvas content */}
  </Canvas>
</div>
```

**Drop Overlay Component**:
- Semi-transparent overlay with upload icon
- Shows supported formats
- Animates on drag enter/leave
- Material-UI styled

#### B. File Handling Logic

```typescript
const handleDrop = async (e: DragEvent) => {
  e.preventDefault();
  const files = Array.from(e.dataTransfer.files);

  for (const file of files) {
    // Check if file is supported
    const ext = file.name.split('.').pop()?.toLowerCase();
    if (!isSupportedFormat(ext)) {
      showError(`Unsupported format: ${ext}`);
      continue;
    }

    // Upload file content
    await uploadFileContent(file);
  }
};
```

### Backend Implementation

#### A. New Upload Endpoint (`/api/file-browser/upload`)

```python
@file_browser.route("/upload", methods=["POST"])
def upload_file():
    """
    Upload file content for loading into ZnDraw.

    Form Data:
    - file: File to upload (multipart/form-data)
    - room: Optional room name
    - start, stop, step: Optional slice parameters
    - make_default: Optional boolean
    """
    if 'file' not in request.files:
        return jsonify({"error": "No file provided"}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({"error": "No file selected"}), 400

    # Validate file type
    ext = Path(file.filename).suffix.lstrip('.').lower()
    if ext not in FORMAT_BACKENDS:
        return jsonify({"error": f"Unsupported format: {ext}"}), 400

    # Save to temporary location
    temp_dir = current_app.config.get("UPLOAD_TEMP_DIR", "/tmp/zndraw_uploads")
    os.makedirs(temp_dir, exist_ok=True)

    # Generate unique filename to avoid collisions
    unique_id = str(uuid.uuid4())
    temp_filename = f"{unique_id}_{secure_filename(file.filename)}"
    temp_path = os.path.join(temp_dir, temp_filename)

    file.save(temp_path)

    # Get or generate room name
    room = request.form.get("room")
    if not room:
        room = generate_room_name(file.filename, redis_client)

    # Queue Celery task with temp file path
    task = read_file.delay(
        file=temp_path,
        room=room,
        server_url=current_app.config.get("SERVER_URL"),
        start=request.form.get("start", type=int),
        stop=request.form.get("stop", type=int),
        step=request.form.get("step", type=int),
        make_default=request.form.get("make_default", type=bool, default=False),
        cleanup_after=True,  # NEW: Delete temp file after loading
    )

    return jsonify({
        "status": "queued",
        "room": room,
        "task_id": task.id,
        "message": "File uploaded and loading queued"
    })
```

#### B. Enhanced Celery Task

```python
@shared_task
def read_file(
    file: str,
    room: str,
    server_url: str = "http://localhost:5000",
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    make_default: bool = False,
    batch_size: int = 10,
    root_path: str | None = None,
    cleanup_after: bool = False,  # NEW parameter
) -> None:
    try:
        # ... existing file loading logic ...

        # Store metadata about upload source
        if cleanup_after:
            metadata_manager.update({
                "upload_source": "browser_upload",
                "original_filename": Path(file).name
            })

    finally:
        # Cleanup temporary file if requested
        if cleanup_after and Path(file).exists():
            try:
                os.remove(file)
                vis.log(f"Cleaned up temporary file")
            except Exception as e:
                log.error(f"Failed to cleanup temp file: {e}")
```

### Frontend API Client

```typescript
// myapi/client.ts
export interface UploadFileRequest {
  file: File;
  room?: string;
  start?: number;
  stop?: number;
  step?: number;
  make_default?: boolean;
}

export const uploadFile = async (request: UploadFileRequest): Promise<{
  status: string;
  room: string;
  task_id: string;
  message: string;
}> => {
  const formData = new FormData();
  formData.append('file', request.file);

  if (request.room) formData.append('room', request.room);
  if (request.start !== undefined) formData.append('start', request.start.toString());
  if (request.stop !== undefined) formData.append('stop', request.stop.toString());
  if (request.step !== undefined) formData.append('step', request.step.toString());
  if (request.make_default) formData.append('make_default', 'true');

  const { data } = await apiClient.post('/api/file-browser/upload', formData, {
    headers: { 'Content-Type': 'multipart/form-data' }
  });

  return data;
};
```

---

## Feature 2: Copy/Paste Support

### Use Cases

1. **File Copy from File Manager**: User copies file, pastes in canvas
2. **Text Copy**: User copies file path as text, pastes in canvas
3. **URL Copy**: User copies URL to remote file (future enhancement)

### Implementation

#### A. Clipboard Event Handler

```typescript
// Add to Canvas or App component
useEffect(() => {
  const handlePaste = async (e: ClipboardEvent) => {
    e.preventDefault();

    // Handle file paste
    const files = Array.from(e.clipboardData?.files || []);
    if (files.length > 0) {
      for (const file of files) {
        await handleFileUpload(file);
      }
      return;
    }

    // Handle text paste (file paths)
    const text = e.clipboardData?.getData('text');
    if (text) {
      await handleTextPaste(text);
    }
  };

  window.addEventListener('paste', handlePaste);
  return () => window.removeEventListener('paste', handlePaste);
}, []);
```

#### B. Text Path Handler

```typescript
const handleTextPaste = async (text: string) => {
  // Check if text looks like a file path
  const possiblePaths = text.split('\n').filter(line =>
    line.trim() && (
      line.includes('/') ||
      line.includes('\\') ||
      line.match(/\.[a-z0-9]+$/i)
    )
  );

  if (possiblePaths.length === 0) return;

  for (const path of possiblePaths) {
    const trimmedPath = path.trim();

    // Try to load via file browser (requires server access)
    try {
      await loadFile({ path: trimmedPath });
    } catch (error) {
      showError(`Cannot load path: ${trimmedPath}`);
    }
  }
};
```

#### C. Copy Support (Ctrl+C)

```typescript
// Copy current scene/selection
const handleCopy = async (e: ClipboardEvent) => {
  if (!roomId) return;

  const selection = useAppStore.getState().selection;

  if (selection && selection.length > 0) {
    // Copy selected particles as XYZ text
    const xyzText = await exportSelectionAsXYZ(roomId, selection);
    e.clipboardData?.setData('text/plain', xyzText);
    e.preventDefault();
  }
};
```

---

## Feature 3: File Download

### Use Cases

1. **Download Current Frame**: Export single frame as XYZ, PDB, etc.
2. **Download Trajectory**: Export selected frames or full trajectory
3. **Download Selection**: Export only selected particles
4. **Download Geometry**: Export geometry objects (curves, planes, etc.)

### Implementation

#### A. Frontend Download UI

```typescript
// Add to sidebar or context menu
<Menu>
  <MenuItem onClick={() => handleDownload('current', 'xyz')}>
    Download Current Frame (XYZ)
  </MenuItem>
  <MenuItem onClick={() => handleDownload('current', 'pdb')}>
    Download Current Frame (PDB)
  </MenuItem>
  <MenuItem onClick={() => handleDownload('trajectory', 'xyz')}>
    Download Full Trajectory (XYZ)
  </MenuItem>
  <MenuItem onClick={() => handleDownload('selection', 'xyz')}>
    Download Selection (XYZ)
  </MenuItem>
</Menu>
```

#### B. Backend Download Endpoint

```python
@main.route("/api/rooms/<string:room_id>/download", methods=["GET"])
def download_frames(room_id: str):
    """
    Download frames in specified format.

    Query Parameters:
    - format: File format (xyz, pdb, cif, etc.)
    - frames: Comma-separated frame indices or 'all'
    - selection: Optional comma-separated particle indices
    - filename: Optional custom filename
    """
    from zndraw.io import ZarrIO
    import ase.io
    from io import BytesIO

    format_type = request.args.get('format', 'xyz')
    frames_param = request.args.get('frames', 'current')
    selection_param = request.args.get('selection')
    custom_filename = request.args.get('filename')

    # Validate format
    if format_type not in FORMAT_BACKENDS:
        return jsonify({"error": f"Unsupported format: {format_type}"}), 400

    # Get frames
    io = ZarrIO(room_id=room_id, redis_client=redis_client)

    if frames_param == 'all':
        atoms_list = list(io)
    elif frames_param == 'current':
        current_frame = redis_client.get(f"room:{room_id}:frame")
        atoms_list = [io[int(current_frame)]]
    else:
        frame_indices = [int(i) for i in frames_param.split(',')]
        atoms_list = [io[i] for i in frame_indices]

    # Apply selection filter
    if selection_param:
        selection = [int(i) for i in selection_param.split(',')]
        atoms_list = [atoms[selection] for atoms in atoms_list]

    # Write to buffer
    buffer = BytesIO()

    if len(atoms_list) == 1:
        ase.io.write(buffer, atoms_list[0], format=format_type)
    else:
        ase.io.write(buffer, atoms_list, format=format_type)

    buffer.seek(0)

    # Generate filename
    if not custom_filename:
        if len(atoms_list) == 1:
            custom_filename = f"{room_id}_frame.{format_type}"
        else:
            custom_filename = f"{room_id}_trajectory.{format_type}"

    return send_file(
        buffer,
        mimetype='application/octet-stream',
        as_attachment=True,
        download_name=custom_filename
    )
```

#### C. Frontend Download Function

```typescript
// myapi/client.ts
export const downloadFrames = async (
  roomId: string,
  format: string,
  frames: 'current' | 'all' | number[],
  selection?: number[],
  filename?: string
): Promise<void> => {
  const params = new URLSearchParams();
  params.append('format', format);

  if (frames === 'current' || frames === 'all') {
    params.append('frames', frames);
  } else {
    params.append('frames', frames.join(','));
  }

  if (selection) {
    params.append('selection', selection.join(','));
  }

  if (filename) {
    params.append('filename', filename);
  }

  // Trigger browser download
  const url = `/api/rooms/${roomId}/download?${params.toString()}`;
  window.location.href = url;
};
```

---

## UI/UX Enhancements

### 1. Upload Progress Indicator

```typescript
// UploadProgress component
<Dialog open={uploadInProgress}>
  <DialogContent>
    <CircularProgress />
    <Typography>Uploading {currentFile}...</Typography>
    <LinearProgress variant="determinate" value={progress} />
  </DialogContent>
</Dialog>
```

### 2. Format Detection Feedback

```typescript
// Show format info on hover during drag
<DropOverlay>
  <Typography variant="h4">Drop file to upload</Typography>
  {draggedFormat && (
    <Chip
      label={`${draggedFormat} - Supported by ${getBackend(draggedFormat)}`}
      color="success"
    />
  )}
</DropOverlay>
```

### 3. Keyboard Shortcuts

```typescript
// Add to KeyboardShortcutsHandler.tsx
const shortcuts = {
  'Ctrl+V': handlePaste,
  'Ctrl+C': handleCopy,
  'Ctrl+S': () => handleDownload('current', 'xyz'),
  'Ctrl+Shift+S': () => handleDownload('all', 'xyz'),
};
```

---

## Security Considerations

1. **File Size Limits**
   ```python
   MAX_UPLOAD_SIZE = 500 * 1024 * 1024  # 500 MB
   app.config['MAX_CONTENT_LENGTH'] = MAX_UPLOAD_SIZE
   ```

2. **File Type Validation**
   - Validate extension against FORMAT_BACKENDS
   - Check magic bytes for file type verification
   - Reject executable files

3. **Temp File Cleanup**
   - Automatic cleanup after task completion
   - Scheduled cleanup of old temp files
   - Disk quota monitoring

4. **Path Traversal Protection**
   - Use `secure_filename()` from werkzeug
   - Validate temp directory access
   - Sandbox temp file operations

---

## Implementation Phases

### Phase 1: Basic Drag & Drop (Week 1)
- ✅ Canvas drop zone UI
- ✅ File upload endpoint
- ✅ Enhanced Celery task with cleanup
- ✅ Format validation
- ✅ Basic progress indication

### Phase 2: Copy/Paste (Week 2)
- ✅ Paste event handling
- ✅ File paste support
- ✅ Text path paste support
- ✅ Copy selection as XYZ

### Phase 3: Download (Week 3)
- ✅ Download endpoint
- ✅ Format selection UI
- ✅ Frame range selection
- ✅ Selection filter support

### Phase 4: Polish & Advanced Features (Week 4)
- ✅ Upload progress tracking
- ✅ Duplicate detection for uploads
- ✅ Keyboard shortcuts
- ✅ Error handling & user feedback
- ✅ Documentation

---

## Testing Strategy

### Unit Tests
- File format validation
- Path security checks
- Temp file cleanup

### Integration Tests
- Upload → Celery → Loading flow
- Duplicate detection
- Progress tracking

### E2E Tests
- Drag & drop XYZ file
- Copy/paste file
- Download current frame
- Download trajectory

---

## Configuration

```python
# Add to server config
class Config:
    # Drag & Drop / Upload
    UPLOAD_TEMP_DIR = os.getenv("ZNDRAW_UPLOAD_TEMP", "/tmp/zndraw_uploads")
    MAX_UPLOAD_SIZE = int(os.getenv("ZNDRAW_MAX_UPLOAD_MB", "500")) * 1024 * 1024

    # Temp file cleanup
    TEMP_FILE_RETENTION_HOURS = int(os.getenv("ZNDRAW_TEMP_RETENTION_HOURS", "24"))
```

---

## Success Metrics

- ✅ User can drag XYZ file onto canvas → file loads
- ✅ User can paste file → file loads
- ✅ User can download current frame in multiple formats
- ✅ User can download trajectory
- ✅ Upload progress is visible
- ✅ Errors are clearly communicated
- ✅ Temp files are cleaned up automatically

---

## Future Enhancements

1. **Remote URL Support**: Paste URL → download → load
2. **Batch Upload**: Multiple files at once with queue management
3. **Cloud Storage**: Direct integration with S3, Dropbox, etc.
4. **Format Conversion**: Auto-convert unsupported formats
5. **Drag & Drop Reordering**: Drag frames to reorder trajectory
6. **Export Presets**: Save download format preferences
