# Array Editor for Geometry Forms

This document explains the new array editing functionality for geometry forms.

## Overview

The array editor allows users to:
1. **Edit arrays inline** using a MUI-X DataGrid table interface
2. **Snapshot dynamic values** (like `arrays.positions`) to static arrays for persistence
3. **Switch between modes**: dropdown → array → dropdown
4. **Handle single-value and per-instance data** intelligently

## Components

### 1. `ArrayEditorDialog.tsx`
Modal dialog with MUI-X DataGrid for editing `number[]` or `number[][]` values.

**Features**:
- Add/delete rows with buttons
- Inline cell editing with validation
- Paste from clipboard (CSV, TSV, JSON)
- Single-value mode for scalar properties (radius, scale)
- Switch back to dropdown selection
- Real-time validation with error messages

### 2. `ArrayFieldToolbar.tsx`
Toolbar with action buttons displayed next to form fields.

**Buttons**:
- **Snapshot** (camera icon): Convert dynamic reference to static array
- **Edit** (table icon): Open table editor for arrays/numbers
- **Clear** (X icon): Remove static value and return to dropdown

### 3. `arrayEditor.ts` (Utils)
Utility functions for array handling.

**Key Functions**:
- `normalizeToArray()`: Convert any value to 2D array format
- `denormalizeFromArray()`: Convert back to appropriate type
- `validateArrayData()`: Validate dimensions and ranges
- `inferFieldType()`: Auto-detect field type from path
- `parseClipboardData()`: Parse pasted data

### 4. Modified `DynamicEnumRenderer.tsx`
Enhanced to support array editing alongside existing autocomplete.

## Field Types

The system automatically detects field types based on the property name:

| Field Type | Dimensions | Columns | Single Value Support |
|------------|------------|---------|---------------------|
| `position` | 3 | X, Y, Z | No |
| `direction` | 3 | X, Y, Z | No |
| `color` | 3 | R, G, B | No |
| `radius` | 1 | Radius | **Yes** |
| `scale` | 1 | Scale | **Yes** |
| `rotation` | 3 | X, Y, Z | No |

## Single-Value Mode (Radius & Scale)

For fields that support single-value mode:

- **One row** = applies to ALL instances (like `radius: 1.0`)
- **Multiple rows** = per-instance values (like `radius: [1.0, 0.8, 1.2, ...]`)
- Toggle with the chip in the table toolbar

## User Workflows

### Workflow 1: Make Dynamic Value Persistent
1. Field shows `arrays.positions` (dynamic reference)
2. Click **snapshot** button (camera icon)
3. System fetches current frame data
4. Value converted to static array `[[x,y,z], [x,y,z], ...]`
5. Click **edit** button to modify in table

### Workflow 2: Edit Existing Array
1. Field shows `[[1,0,0], [0,1,0]]` (static array)
2. Click **edit** button (table icon)
3. Table editor opens with 2 rows
4. Edit cells, add/delete rows
5. Click **Save**

### Workflow 3: Switch Back to Dropdown
1. In table editor
2. Click **Switch to Dropdown** button
3. Select value from autocomplete (e.g., `arrays.positions`)
4. Static array replaced with dynamic reference

### Workflow 4: Single Value for All Instances
1. Field is `radius` with value `1.0`
2. Click **edit** button
3. Table shows 1 row (applies to all)
4. Click **Add Row** to switch to per-instance mode
5. Each row now represents one instance

## Data Formats

### Input Formats (Normalized to 2D Array)
```typescript
// Single number (for radius/scale)
1.0 → [[1.0]]

// 1D array (multiple instances of scalar)
[1.0, 0.8, 1.2] → [[1.0], [0.8], [1.2]]

// 1D array (single instance of vector)
[1, 0, 0] → [[1, 0, 0]]

// 2D array (multiple instances of vector)
[[1,0,0], [0,1,0]] → [[1,0,0], [0,1,0]]
```

### Output Formats (Denormalized)
```typescript
// Single row with single-value support
[[1.0]] → 1.0

// Single row without single-value support
[[1,0,0]] → [[1,0,0]]

// Multiple rows for 1D field (radius)
[[1.0], [0.8]] → [1.0, 0.8]

// Multiple rows for 3D field (position)
[[1,0,0], [0,1,0]] → [[1,0,0], [0,1,0]]
```

## Clipboard Paste

Supports multiple formats:

**CSV**:
```
1.0, 0.0, 0.0
0.0, 1.0, 0.0
```

**TSV** (Tab-separated):
```
1.0    0.0    0.0
0.0    1.0    0.0
```

**JSON**:
```json
[[1,0,0], [0,1,0]]
```

## Validation

The system validates:
- **Dimensions**: All rows must have correct number of columns
- **Value ranges**: E.g., colors must be 0-1
- **Consistency**: Instance counts should match across fields

Validation errors are shown in the dialog and prevent saving.

## Integration with JSONForms

The array editor integrates seamlessly with the existing JSONForms system:

1. Matches fields with `x-custom-type="dynamic-enum"`
2. Automatically handles both string and array values
3. Preserves data through JSONForms `handleChange`
4. Works with auto-save (250ms debounce)

## Future Enhancements

Potential improvements:
- Drag-and-drop row reordering
- Copy/paste individual rows
- Formula support (e.g., `=row*0.1`)
- Visualization preview (for positions, colors)
- Undo/redo
- Import from file (CSV, JSON)
- Export to file
